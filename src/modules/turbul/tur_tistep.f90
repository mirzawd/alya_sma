!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_tistep()
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_tistep
  ! NAME 
  !    tur_tittim
  ! DESCRIPTION
  !    This routine sets the time step
  ! USES
  ! USED BY
  !    tur_begite
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_turbul
  use mod_outfor, only : outfor
  use mod_communications, only : PAR_MAX
  implicit none

  if(kfl_timco/=2) then   ! not LOCAL time step


     dtinv_tur=dtinv
     if(kfl_stead_tur==1) dtinv_tur = 0.0_rp
     if(kfl_timei_tur==0) dtinv_tur = 0.0_rp 
     !     if( kfl_sgsti_nsi == 0) dtsgs_nsi = 0.0_rp   SGS not ready in turbul
     
     if(kfl_tisch_tur==1) then
        !
        ! Trapezoidal rule: Euler iterations
        !     
        pabdf_tur(1) = 1.0_rp
        pabdf_tur(2) = -1.0_rp        
        if(ittim<=neule_tur) then
           kfl_tiacc_tur=1
        else
           kfl_tiacc_tur=kfl_tiaor_tur
        end if
        if(kfl_tiacc_tur==2) then
!           dtinv_tur = 2.0_rp*dtinv_tur
           pabdf_tur(1) =  2.0_rp  ! from jul 2016 theta is no longer included in dtinv_tur , instead pabdf_tur is exteded to CN also
           pabdf_tur(2) = -2.0_rp
        end if
        nbdfp_tur = 2
        
     else if(kfl_tisch_tur==2) then
        !
        ! BDF scheme: increase integration order at each time step
        !
        if(ittim<=neule_tur) then
           kfl_tiacc_tur=1
        else
           kfl_tiacc_tur=min(kfl_tiaor_tur,ittim)
        end if
        call parbdf(kfl_tiacc_tur,pabdf_tur)
        nbdfp_tur = kfl_tiacc_tur + 1
     end if

  else

     nbdfp_tur = 2       ! This line is needed the 4 below I am not totally sure
     pabdf_tur(1) = 1.0_rp
     pabdf_tur(2) = -1.0_rp
     if( kfl_tiacc_tur > 1 ) call runend('LOCAL TIME STEP WITH SECOND ORDER METHODS NEEDS THINKING')
     
  end if
  

  routp(1)=dtcri_tur

  ! minimum time step when using local time step
  routp(2)=dtcri_tur*max(saflo_tur, safet_tur)
  !
  ! Look for maximum (time step) over subdomains
  !
  if(kfl_timco==2.and.kfl_paral>=0)  &                      
       call PAR_MAX(dtmax_tur)
  routp(3)=dtmax_tur
  ioutp(1)=kfl_timei_tur
  ioutp(2)=kfl_stead_tur
  call outfor(8_ip,lun_outpu,' ')

end subroutine tur_tistep
