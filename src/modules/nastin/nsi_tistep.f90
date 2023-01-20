!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_tistep
!-----------------------------------------------------------------------
!****f* Nastin/nsi_tistep
! NAME 
!    nsi_tistep
! DESCRIPTION
!    This routine sets the time step. 
! USES
! USED BY
!    nsi_begite
!***
!-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_nastin
  use mod_outfor, only : outfor
  use mod_communications,        only : PAR_MAX
  implicit none

  if(kfl_timco/=2) then

     dtinv_nsi = dtinv
     if( kfl_timei_nsi == 0) dtinv_nsi = 0.0_rp 
     if( kfl_stead_nsi == 1) dtinv_nsi = 0.0_rp
     if( kfl_sgsti_nsi == 0) dtsgs_nsi = 0.0_rp

     if(kfl_tisch_nsi==1) then
        !
        ! Trapezoidal rule: Euler iterations
        !
        pabdf_nsi(1) = 1.0_rp
        pabdf_nsi(2) = -1.0_rp
        if(ittim<=neule_nsi) then
           kfl_tiacc_nsi = 1
        else
           kfl_tiacc_nsi = kfl_tiaor_nsi
        end if
        dtsgs_nsi = dtinv_nsi
        if( kfl_tiacc_nsi == 2 ) then
           pabdf_nsi(1) =  2.0_rp  ! from vers 772 theta is no longer included in dtinv_nsi , instead pabdf_nsi is exteded to CN also
           pabdf_nsi(2) = -2.0_rp
           if( kfl_sgsac_nsi == 2 ) dtsgs_nsi = 2.0_rp*dtsgs_nsi ! vers 772 sgs CN need to be updated for programing consistency
        end if
        nbdfp_nsi = 2

     else if(kfl_tisch_nsi==2) then
        !
        ! BDF scheme: increase integration order at each time step
        !
        dtsgs_nsi = dtinv_nsi
        if(ittim<=neule_nsi) then
           kfl_tiacc_nsi=1
        else
           kfl_tiacc_nsi=min(kfl_tiaor_nsi,ittim)
        end if
        call parbdf(kfl_tiacc_nsi,pabdf_nsi)
        nbdfp_nsi = kfl_tiacc_nsi + 1

     else if(kfl_tisch_nsi==3) then
        !
        ! Explicit Adams-Bashforth schemme
        !
        ! k1*u^{n+1}  = k2*u^{n} + M^{-1}*\Delta_t * [k3*rm^{n} + k4*rm^{n-1} + k5*rm^{n-2}] 

        pabdf_nsi(1) = 1.0_rp
        pabdf_nsi(2) = -1.0_rp
        if     (kfl_tiacc_nsi == 1) then
           pabdf_nsi(3) = 1.0_rp
           pabdf_nsi(4) = 0.0_rp
           pabdf_nsi(5) = 0.0_rp
        else if(kfl_tiacc_nsi == 2) then
           pabdf_nsi(3) = 3.0_rp/2.0_rp
           pabdf_nsi(4) = -1.0_rp/2.0_rp
           pabdf_nsi(5) = 0.0_rp
        else if(kfl_tiacc_nsi == 3) then
           pabdf_nsi(3) = 23.0_rp/12.0_rp
           pabdf_nsi(4) = -16.0_rp/12.0_rp
           pabdf_nsi(5) = 5.0_rp/12.0_rp
        else
           call runend('TISTEP: AB not coded for order > 3')
        end if
        nbdfp_nsi = 2
        dtsgs_nsi = dtinv_nsi
     else if(kfl_tisch_nsi==4) then
        pabdf_nsi(1) = 3.0_rp/2.0_rp
        pabdf_nsi(2) = -4.0_rp/2.0_rp
        pabdf_nsi(3) = 1.0_rp/2.0_rp
        nbdfp_nsi = 3
        dtsgs_nsi = dtinv_nsi
     end if

     !
     ! If dt is prescribed, compute corresponding safety factor
     !
     if(kfl_timco==0) then
        if(dtcri_nsi/=0.0_rp.and.dtinv_nsi/=0.0_rp) then
           safet_nsi=1.0_rp/(dtcri_nsi*dtinv_nsi)
        else
           safet_nsi=1.0_rp
        end if
     end if
 
  else  ! local time step

     nbdfp_nsi = 2       ! This line is needed the 4 below I am not totally sure
     pabdf_nsi(1) = 1.0_rp
     pabdf_nsi(2) = -1.0_rp
     dtsgs_nsi = dtinv_nsi
     if( kfl_tiacc_nsi > 1 ) call runend('LOCAL TIME STEP WITH SECOND ORDER METHODS NEEDS THINKING')

  end if

  routp(1) = dtcri_nsi
  ioutp(1) = kfl_timei_nsi
  ioutp(2) = kfl_stead_nsi

  ! minimum time step when using local time step
  routp(2) = max(saflo_nsi*sqrt(safet_nsi/safma_nsi), safet_nsi)*dtcri_nsi 

  !
  ! Look for maximum  (time step) over subdomains
  !
  if(kfl_timco==2.and.kfl_paral>=0)  &                      
       call PAR_MAX(dtmax_nsi)
  routp(3) = dtmax_nsi

  call outfor(8_ip,lun_outpu,' ')

end subroutine nsi_tistep
