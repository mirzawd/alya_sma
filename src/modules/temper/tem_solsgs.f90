!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_solsgs()
  !-----------------------------------------------------------------------
  !****f* Temper/tem_solsgs
  ! NAME 
  !    tem_solsgs
  ! DESCRIPTION
  !    This routine solves the SGS equation
  ! USES
  ! USED BY
  !    tem_endite
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_domain
  use def_master
  use def_temper
  use mod_memory
  use mod_ADR,     only : PROJECTIONS_AND_SGS_ASSEMBLY ! 4
  use mod_ADR,     only : BUBBLE_ASSEMBLY              ! 5
  use mod_ADR,     only : ADR_initialize_projections
  use mod_ADR,     only : ADR_update_projections
  implicit none
  real(rp)    :: time1,time2
   
  if(  ADR_tem % kfl_time_sgs      /= 0 .or. &
       ADR_tem % kfl_nonlinear_sgs /= 0 .or. &
       ADR_tem % kfl_stabilization  > 0 ) then
     
     !-------------------------------------------------------------------
     !
     ! Update subgrid scale and project residuals
     !
     !-------------------------------------------------------------------

     call cputim(time1)
     resgs_tem(1) = 0.0_rp
     resgs_tem(2) = 0.0_rp
     itsgm_tem    = 0
     !
     ! Residual projections
     !
     if( INOTMASTER ) then
        !
        ! Initialize
        !
        call ADR_initialize_projections(ADR_tem)
        !
        ! Update SGS and projections
        !              
        call tem_elmope_new(PROJECTIONS_AND_SGS_ASSEMBLY)
        !
        ! Residual projections
        ! 
        call ADR_update_projections(ADR_tem,vmass)
     end if
 
     call cputim(time2)    
     !
     ! Output SGS convergence
     !
     call tem_cvgsgs()

  else if( ADR_tem % kfl_time_bubble /= 0 ) then

     !-------------------------------------------------------------------
     !
     ! Update bubble
     !
     !-------------------------------------------------------------------

     if( INOTMASTER ) call tem_elmope_new(BUBBLE_ASSEMBLY)

  end if

end subroutine tem_solsgs
