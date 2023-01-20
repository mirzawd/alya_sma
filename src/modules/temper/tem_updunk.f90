!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Temper
!> @{
!> @file    tem_updunk.f90
!> @author  houzeaux
!> @date    2018-12-28
!> @brief   This routine performs several types of updates for the
!>          thermal variable (T or h)
!> @details Solution updates
!> @}
!-----------------------------------------------------------------------

subroutine tem_updunk(itask)

  use def_parame
  use def_master
  use def_domain
  use def_temper
  use mod_ADR,    only : ADR_end_time_step
  use mod_ADR,    only : ADR_begin_inner_iteration
  use mod_ADR,    only : ADR_end_inner_iteration
  use mod_ADR,    only : ADR_begin_time_step
  use mod_ADR,    only : ADR_after_restart
  use mod_arrays, only : arrays_number
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: itime,ipoin
  integer(ip)             :: nprev_tem,icomp
  real(rp)                :: rela1_tem

  nprev_tem = min(3_ip,ncomp_tem)
  
  select case (itask)

  case ( ITASK_INIUNK )
     !
     ! (:,all) <= (:,1): Initial solution
     !
     do icomp = 1,ncomp_tem
        do ipoin = 1,npoin
           therm(ipoin,icomp) = therm(ipoin,1)
        end do
     end do
     
  case ( ITASK_BEGSTE )
     !
     ! (:,2) <= (:,1): Initial guess for outer iterations
     !
     if (INOTMASTER) then ! avoid null ptr with nvfortran -g on power9
       call ADR_begin_time_step(ADR_tem,therm)
     endif
     if( kfl_regim_tem == 4 ) then
        do ipoin=1,nunkn_tem
           tempe(ipoin,2) = tempe(ipoin,1)
        end do
     endif
     do ipoin = 1,npoin
        therm(ipoin,2) = therm(ipoin,1)
     end do

  case ( ITASK_BEGITE )
     !
     ! UNKNO <= (:,1)
     !
     do ipoin=1,nunkn_tem
        unkno(ipoin) = therm(ipoin,1)
     end do

  case ( ITASK_INNITE )
     !
     ! (:,1) <= UNKNO: Inside Runge-Kutta 
     !
     if(postp(1) % npp_stepi(arrays_number('RESID'),0)>0) then
        do ipoin=1,nunkn_tem
           teold_tem(ipoin) = therm(ipoin,1)
        end do
     end if
     if(relax_tem==1.0_rp) then
        do ipoin=1,nunkn_tem
           therm(ipoin,1) = unkno(ipoin)
        end do
     else
        rela1_tem=1.0_rp-relax_tem
        do ipoin=1,nunkn_tem
           therm(ipoin,1) = relax_tem*unkno(ipoin)+rela1_tem*therm(ipoin,1)
        end do
     end if
     
  case ( ITASK_ENDINN )
     !
     ! (:,1) <= UNKNO
     !
     if(postp(1) % npp_stepi(arrays_number('RESID'),0)>0) then
        do ipoin=1,nunkn_tem
           teold_tem(ipoin) = therm(ipoin,1)
        end do
     end if
     if(relax_tem==1.0_rp) then
        do ipoin=1,nunkn_tem
           therm(ipoin,1) = unkno(ipoin)
        end do
     else
        rela1_tem=1.0_rp-relax_tem
        do ipoin=1,nunkn_tem
           therm(ipoin,1) = relax_tem*unkno(ipoin)+rela1_tem*therm(ipoin,1)
        end do
     end if
     
  case ( ITASK_ENDITE )
     !
     ! (:,2) <= (:,1): End of inner iteration
     !        
     if (INOTMASTER) then ! avoid null ptr with nvfortran -g on power9
       call ADR_end_inner_iteration(ADR_tem,therm)
     end if
     call tem_temperature_from_enthalpy()
     if( kfl_regim_tem == 4_ip .and. kfl_lookg_tem <= 0 ) then        
        do ipoin = 1,npoin
           tempe(ipoin,2) = tempe(ipoin,1)
        end do
     end if
     
  case ( ITASK_ENDSTE )
     !
     ! (:,3) <= (:,1): End of time step
     ! (:,4) <= (:,3)
     ! (:,5) <= (:,4)
     ! ...
     !        
     if (INOTMASTER) then ! avoid null ptr with nvfortran -g on power9
       call ADR_end_time_step(ADR_tem,therm)
     end if
     call tem_temperature_from_enthalpy()
     
     if( kfl_regim_tem == 4_ip .and. kfl_lookg_tem <= 0 ) then
        ! 
        ! High-order temporal schemes 
        !
        if( kfl_tisch_tem == 2 ) then
           do itime = 2+kfl_tiaor_tem,4,-1
              do ipoin=1,npoin
                 tempe(ipoin,itime) = tempe(ipoin,itime-1)
              end do
           end do
        end if
        do ipoin = 1,npoin
           tempe(ipoin,nprev_tem) = tempe(ipoin,1)
        end do
     end if
             
  end select


end subroutine tem_updunk
