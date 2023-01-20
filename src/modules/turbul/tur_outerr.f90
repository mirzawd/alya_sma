!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_outerr
  !------------------------------------------------------------------------
  !****f* Temper/tur_outerr
  ! NAME 
  !    tur_outerr
  ! DESCRIPTION
  !    This routine checks if there are errros and warnings
  ! USES
  ! USED BY
  !    tur_turnon
  !***
  !------------------------------------------------------------------------
  use      def_master
  use      def_turbul
  use mod_outfor, only : outfor
  use mod_output_postprocess, only : output_postprocess_check_variable_postprocess
  use mod_arrays,             only : arrays_number
  implicit none
  integer(ip)   :: ierro(1)=0,iwarn(1)=0
  !
  ! Check the transient evolution
  !
  if(kfl_timei/=0) then
     if(kfl_timei_tur == 0) then
        iwarn=iwarn+1
        call outfor(2_ip,momod(modul)%lun_outpu,&
             'STEADY EQUATION IN A TRANSIENT CALCULATION')
     end if
  end if
  if(solve(1)%kfl_algso==-2.and.kfl_timei_tur==0) then
     ierro=ierro+1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'EXPLICIT SOLVER REQUIRES TIME TERM TO BE ON')     
  end if
  !
  ! Check postprocess requests
  !
  if(TUR_SPALART_ALLMARAS) then
     if( output_postprocess_check_variable_postprocess(arrays_number('KEY  ')) ) then
        iwarn=iwarn+1
        call outfor(2_ip,momod(modul)%lun_outpu,&
             'CANNOT POSTPROCESS KEY')
        postp(1) % npp_stepi(arrays_number('KEY  '),:)=0
        postp(1) % pos_times(:,arrays_number('KEY  '),:)=0.0_rp
     end if
     if( output_postprocess_check_variable_postprocess(arrays_number('EPSIL')) ) then
        iwarn=iwarn+1
        call outfor(2_ip,momod(modul)%lun_outpu,&
             'CANNOT POSTPROCESS EPSILON')
        postp(1) % npp_stepi(arrays_number('EPSIL'),:)=0
        postp(1) % pos_times(:,arrays_number('EPSIL'),:)=0.0_rp
     end if
     if( output_postprocess_check_variable_postprocess(arrays_number('OMEGA')) ) then
        iwarn=iwarn+1
        call outfor(2_ip,momod(modul)%lun_outpu,&
             'CANNOT POSTPROCESS OMEGA')
        postp(1) % npp_stepi(arrays_number('OMEGA'),:)=0
        postp(1) % pos_times(:,arrays_number('OMEGA'),:)=0.0_rp
     end if
  end if
  if(nturb_tur==1.and. output_postprocess_check_variable_postprocess(arrays_number('RESI2')) ) then
     iwarn=iwarn+1
     call outfor(2_ip,momod(modul)%lun_outpu,&
          'CANNOT POSTPROCESS SECOND TURBULENCE VARIABLE DIFFERENCE')
     postp(1) % npp_stepi(arrays_number('RESI2'),:)=0
     postp(1) % pos_times(:,arrays_number('RESI2'),:)=0.0_rp     
  end if
  if(.not.TUR_SPALART_ALLMARAS) then
     if( output_postprocess_check_variable_postprocess(arrays_number('NUTIL')) ) then
        iwarn=iwarn+1
        call outfor(2_ip,momod(modul)%lun_outpu,&
             'CANNOT POSTPROCESS NU TILDE')
        postp(1) % npp_stepi(arrays_number('NUTIL'),:)=0
        postp(1) % pos_times(:,arrays_number('NUTIL'),:)=0.0_rp
     end if
  end if
  if(kfl_infl2_tur==3.and.(kfl_fixn8_tur==1.or.kfl_fixn6_tur==1)) then
     if(nutnu_tur<0.0_rp) then
        ierro=ierro+1
        call outfor(1_ip,momod(modul)%lun_outpu,&
             'VISCOSITY RATIO NUT/NU IS NEEDED AND MUST BE POSITIVE')
     end if
  end if
  if(kfl_infl1_tur==1.and.(kfl_fixn8_tur==1.or.kfl_fixn6_tur==1)) then 
     if(turin_tur<0.0_rp .and. .not.TUR_SPALART_ALLMARAS) then
        ierro=ierro+1
        call outfor(1_ip,momod(modul)%lun_outpu,&
             'TURBULENCE INTENSITY IS NEEDED AND MUST BE POSITIVE')
     end if
  end if
  if(kfl_algor_tur==2.and.kfl_assem_tur==2) then
     ierro=ierro+1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'CELL ASSEMBLY IMPOSSIBLE WITH COUPLING')     
  end if
  !if(.not.TUR_K_EPS_STD .and. kfl_ortho_tur>=1) then
  !   ierro=ierro+1
  !   call outfor(1_ip,momod(modul)%lun_outpu,&
  !        'OSS AND FULLOSS ARE ONLY WORKING CORRECTLY FOR TUR_K_EPS_STD')     
  !end if

  if ( kfl_sgsti_tur /= 0 ) then
     ierro=ierro+1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'kfl_sgsti_tur /= 0 not ready in turbul - guillaume has added it in reanut but nothing else')     
  end if

  if ( kfl_sgsno_tur /= 0 ) then
     ierro=ierro+1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'kfl_sgsno_tur /= 0 not ready in turbul - guillaume has added it in reanut but nothing else')     
  end if

  if ( kfl_timco == 2 .and. kfl_tiacc_tur == 2 ) then
     ! if eliminated check Lumped mass evolution matrix in tur_elmmsu
     ierro=ierro+1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'Local time step does not make sense with second order time accuracy')     
  end if


  !----------------------------------------------------------------------
  !
  ! ERROR MESSAGE
  !
  !----------------------------------------------------------------------

  call errors(3_ip,ierro,iwarn,' ')

end subroutine tur_outerr
