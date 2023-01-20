!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_begste.f90
!> @author  Solidz Team
!> @date
!> @brief   This routine prepares for a new time step of Solidz
!> @details
!>
!> @}
!-----------------------------------------------------------------------

subroutine sld_begste()

  use def_kintyp,           only : ip
  use def_master,           only : ITASK_BEGSTE, ITASK_ENDSTE
  use def_domain,           only : kfl_elcoh
  use def_solidz,           only : kfl_stead_sld, kfl_rigid_sld
  use def_solidz,           only : kfl_xfeme_sld
  use mod_sld_energy,       only : sld_updene
  use mod_sld_rbo,          only : sld_rbo_updunk, sld_rbo_updbcs
#if defined COMMDOM && COMMDOM == 2
  use def_solidz,           only : kfl_conta_sld
  use mod_sld_pdn_contact_plepp, only : sld_pdn_contact_begste
#endif
  
  implicit none

  if( kfl_stead_sld /= 1_ip ) then

     if( kfl_rigid_sld == 0_ip ) then

        !-------------------------------------------------------------------
        !
        ! Deformable body
        !
        !-------------------------------------------------------------------
        !
        ! Initial guess for the displacement: d(n,0,*) <-- d(n-1,*,*).
        !
        call sld_updunk(ITASK_BEGSTE)  ! Update unknowns (:,ITER_AUX) <= (:,TIME_N)
        call sld_updene(ITASK_BEGSTE)  ! Update energies
        !
        ! Update boundary conditions
        !
        call sld_updbcs(ITASK_BEGSTE)
        !
        ! Updates XFEM/CZM
        !
        if( kfl_elcoh /= 0_ip .or. kfl_xfeme_sld /= 0_ip ) then
           call sld_updcra()           ! Update crack path propagation
           call sld_enrich(1_ip)       ! Enrichement: define enriched nodes
           call sld_uptcoh()           ! Effective traction for cohesive elements
        end if
        !
        ! Coupling: FSI (Use the prediction for the zonal coupling (if any))
        !
        call sld_coupre()
        !
        ! Coupling: PDN contact
        !
#if defined COMMDOM && COMMDOM == 2        
        if( kfl_conta_sld /= 0_ip ) call sld_pdn_contact_begste()
#endif

     else if( kfl_rigid_sld == 1_ip ) then

        !-------------------------------------------------------------------
        !
        ! Rigid body
        !
        !-------------------------------------------------------------------
        !
        ! Initial guess for RB unknowns:
        !
        call sld_rbo_updunk(ITASK_BEGSTE) ! Update unknowns (:,ITER_AUX) <= (:,TIME_N)
        !
        ! Update boundary conditions
        !
        call sld_rbo_updbcs(ITASK_BEGSTE)
        !
        ! Coupling: PDN contact
        !
#if defined COMMDOM && COMMDOM==2   
        if( kfl_conta_sld /= 0_ip ) call sld_pdn_contact_begste()
#endif
        
     end if

  end if

end subroutine sld_begste

