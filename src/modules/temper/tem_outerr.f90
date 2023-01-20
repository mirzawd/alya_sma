!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_outerr()
  !------------------------------------------------------------------------
  !****f* Temper/tem_outerr
  ! NAME 
  !    tem_outerr
  ! DESCRIPTION
  !    This routine checks if there are errros and warnings
  ! USES
  ! USED BY
  !    tem_turnon
  !***
  !------------------------------------------------------------------------
  use def_master
  use def_domain
  use def_temper
  use def_kermod
  use mod_outfor,             only : outfor
  use mod_arrays,             only : arrays_number
  use mod_output_postprocess, only : output_postprocess_check_variable_postprocess
  use mod_output_postprocess, only : output_postprocess_cancel_variable_postprocess
  implicit none
  integer(ip)   :: ierro=0,iwarn=0
  !
  ! Properties
  !
  if( kfl_prope == 0 ) then
     ierro = ierro+1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'WHEN USING TEMPER PROPERTIES SHOUD BE DECALRED IN KERMOD')
  end if
  !
  ! Sources
  !
  if( kfl_sourc_tem == SOURCE_TERM_FIELD .and. INOTMASTER ) then
     if( size(xfiel) < kfl_sonum_tem ) then
        ierro = ierro + 1
        call outfor(1_ip,momod(modul)%lun_outpu,&
             'WRONG FIELD FOR HEAT SOURCE')
     else if( .not. associated(xfiel(kfl_sonum_tem) % a)) then
        ierro = ierro + 1
        call outfor(1_ip,momod(modul)%lun_outpu,&
             'WRONG FIELD FOR HEAT SOURCE')
     else if( kfl_field(2,kfl_sonum_tem) /= NELEM_TYPE ) then
        ierro = ierro + 1
        call outfor(1_ip,momod(modul)%lun_outpu,&
             'WRONG FIELD FOR HEAT SOURCE: FIELD MUST BE OF ELEMENT TYPE')
     end if
  end if
  !
  ! Check the transient evolution
  !
  if( kfl_timei /= 0 ) then
     if(kfl_timei_tem == 0) then
        iwarn=iwarn+1
        call outfor(2_ip,momod(modul)%lun_outpu,&
             'STEADY TEMPERATURE IN A TRANSIENT CALCULATION')
     end if
  end if
  !
  ! Boussinesq without temperature
  !
  if(kfl_advec_tem==1.and.kfl_inter_tem==0.and.kfl_modul(1)==0.and.kfl_modul(6)==0.and.kfl_vefun==0) then
     ierro=ierro+1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'CONVECTION IN TEMPER IS IMPOSSIBLE IF NO OTHER MODULE SOLVE FOR THE VELOCITY')
  end if
  !
  ! Time integration scheme and accuracy
  !
  if(kfl_timei_tem==1.and.kfl_tisch_tem==1.and.kfl_tiacc_tem>2) then
     ierro=ierro+1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'WRONG TIME INTEGRATION ORDER USING TRAPEZOIDAL RULE')
  end if
  !
  ! Time tracking and integration scheme
  !
!!$  if(kfl_timei_tem==1.and.kfl_sgsti_tem/=0.and.kfl_tisch_tem==2) then
!!$     ierro=ierro+1
!!$     call outfor(1_ip,momod(modul)%lun_outpu,&
!!$          'CANNOT TRACK THE SUBGRID SCALES WITH A BDF SCHEME')
!!$  end if
  !
  ! Time tracking of the subscales
  !
  if(kfl_timei_tem==0.and.kfl_sgsti_tem/=0) then
     iwarn=iwarn+1
     kfl_sgsti_tem=0
     call outfor(2_ip,momod(modul)%lun_outpu,&
          'CANNOT TRACK THE SUBGRID SCALES IN TIME FOR STATIONARY PROBLEM')
  end if
  !
  ! Exact solution
  !
  if(kfl_sourc_tem/=0.and.kfl_exacs_tem/=0) then
     iwarn=iwarn+1
     kfl_sourc_tem=0
     call outfor(2_ip,momod(modul)%lun_outpu,&
          'SOURCE TERM WAS AUTOMATICALLY SET TO ZERO TO SOLVE AN EXACT SOLUTION')     
  end if
  !
  ! Turbulence without solving TURBUL
  ! 
!  if(kfl_cotur_tem==1.and.kfl_modul(4)==0.and.kfl_inter_tem==0) then
!     ierro=ierro+1
!     call outfor(1_ip,momod(modul)%lun_outpu,&
!          'TURBULENCE COUPLING IS IMPOSSIBLE IF TURBUL MODULE IS NOT SOLVED')
!  end if
  !
  ! Turbulent flow without turbulent Prandtl numebr
  !
  if(prtur_tem<=0.0_rp.and. turmu_ker % kfl_exist /= 0_ip ) then
     ierro=ierro+1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'PRANDTL NUMBER CANNOT BE NEGATIVE NOR ZERO WHEN COUPLING WITH A TURBULENCE MODEL')
  end if

  !----------------------------------------------------------------------
  !
  ! Postprocess
  !
  !----------------------------------------------------------------------
  !
  ! Orthogonal projection
  !
  if( output_postprocess_check_variable_postprocess(arrays_number('PROJ1')) .and. kfl_ortho_tem == 0 ) then
     call output_postprocess_cancel_variable_postprocess(arrays_number('PROJ1'))
     iwarn = iwarn + 1
     call outfor(2_ip,momod(modul)%lun_outpu,'CANNOT POSTPROCESS ORTHOGONAL PROJECTION')
  end if
  !----------------------------------------------------------------------
  !
  ! ERROR MESSAGE
  !
  !----------------------------------------------------------------------
  call errors(3_ip,ierro,iwarn,' unknown? ')

end subroutine tem_outerr
