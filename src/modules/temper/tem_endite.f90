!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_endite(itask)
!-----------------------------------------------------------------------
!****f* Temper/tem_endite
! NAME 
!    tem_endite
! DESCRIPTION
!    This routine checks convergence and performs updates of the
!    temperature  at:
!    - itask=1 The end of an internal iteration
!    - itask=2 The end of the internal loop iteration
! USES
!    tem_cvgunk
!    tem_updunk
! USED BY
!    tem_doiter
!***
!-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_temper
  use def_kermod,          only : kfl_adj_prob
  use def_kermod,          only : kfl_cos_opt 
  use mod_ADR,             only : ADR_manufactured_error
  use mod_messages,        only : livinf
  use mod_tem_therm_press, only : therm_press_update
  implicit none
  integer(ip), intent(in) :: itask

  select case ( itask )

  case ( ITASK_ENDINN )
     !
     !  Compute convergence residual of the internal iteration (that is,
     !  || T(n,i,j) - T(n,i,j-1)|| / ||T(n,i,j)||) and update unknowns:
     !  T(n,i,j-1) <-- T(n,i,j) 
     !
     call tem_clippi()              ! Cut off undershoots
     call tem_cvgunk(ITASK_ENDINN)  ! Residual:   ||UNKNO(:)-TEMPE(:,1)||
     call tem_updunk(ITASK_ENDINN)  ! Relaxation and update:     TEMPE(:,1)=UNKNO

     if(kfl_regim_tem>=3) then
        ! 
        ! If low Mach Updates thermodynamic pressure 
        !
        call therm_press_update(1_ip)
     end if
     !
     ! solves subgrid scales
     !  
     call tem_solsgs()
     ! 
     ! If low Mach with sgs updates thermodynamic pressure again
     !
     if(kfl_regim_tem>=3.and. kfl_sgsti_tem /= 0) then
        call therm_press_update(1_ip)
     end if
 
  case ( ITASK_ENDITE )
     !
     !  Compute convergence residual of the external iteration (that is,
     !  || T(n,i,*) - T(n,i-1,*)|| / ||T(n,i,*)||) and update unknowns:
     !  T(n,i-1,*) <-- T(n,i,*) 
     !
     call livinf(16_ip,' ',itinn(modul))
     call tem_cvgunk(ITASK_ENDITE)    ! Residual: ||TEMPE(:,2)-TEMPE(:,1)||
     call tem_updunk(ITASK_ENDITE)    ! Update:   TEMPE(:,2) = TEMPE(:,1)

     !!call tem_updhfl()     ! Update:   TFLUX= Low-Mach heat flux contribution

     ! ********************
     !   Test for coupling 
     ! call tem_temexch(2_ip)
     ! call tem_temexch(1_ip)
     ! ********************
     if( kfl_exacs_tem /= 0 .and. kfl_timei_tem == 0 ) then
        call ADR_manufactured_error(ADR_tem,ittim,cutim,therm)
     end if
     
     if (kfl_cos_opt  == 1) call tem_costcal()
     if (kfl_adj_prob == 1) call tem_senscal()

  case ( ITASK_INNITE )
     !
     ! After Runge-Kutta substep:
     ! * Clip under and overshoots.
     ! * THERM(:,1)=UNKNO
     !
     call tem_clippi()
     call tem_updunk(ITASK_INNITE)    ! Relaxation and update:     

  end select

end subroutine tem_endite
