!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_tistep
  !-----------------------------------------------------------------------
  !****f* Levels/lev_tistep
  ! NAME 
  !    lev_tittim
  ! DESCRIPTION
  !    This routine sets the time step
  ! USES
  ! USED BY
  !    lev_begite
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_levels
  use mod_outfor, only : outfor
  implicit none

  if(kfl_timco/=2) then
     dtinv_lev=dtinv
     if(kfl_stead_lev==1) dtinv_lev = 0.0_rp

     if(kfl_timet_lev==1) then
        kfl_tiacc_lev=kfl_tiaor_lev
        call parbdf(kfl_tiacc_lev,pabdf_lev)
     else if(kfl_timet_lev==2) then
        if(kfl_tisch_lev==1) then
           !
           ! Trapezoidal rule: Euler iterations
           !
           if(ittim<=neule_lev) then
              kfl_tiacc_lev=1
           else
              kfl_tiacc_lev=kfl_tiaor_lev
           end if
           if(kfl_tiacc_lev==2) dtinv_lev = 2.0_rp*dtinv_lev
        else
           !
           ! BDF scheme: increase integration order at each time step
           !
           kfl_tiacc_lev=min(kfl_tiaor_lev,ittim)
           call parbdf(kfl_tiacc_lev,pabdf_lev)
        end if
     end if
  else
     dtinv_lev=1.0_rp   ! else it gives uninitialized in L69 lev_elmope
  end if

  routp(1)=dtcri_lev
  routp(2)=0.0_rp
  routp(3)=0.0_rp
  ioutp(1)=1
  ioutp(2)=kfl_stead_lev
  call outfor(8_ip,lun_outpu,' ')

end subroutine lev_tistep
