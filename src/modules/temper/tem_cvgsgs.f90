!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_cvgsgs()
  !-----------------------------------------------------------------------
  !****f* temper/tem_cvgsgs
  ! NAME 
  !    tem_cvgsgs
  ! DESCRIPTION
  !    This routine outputs SGS statistics
  ! USES
  ! USED BY
  !    tem_solsgs 
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_temper
  implicit none
!  integer(ip), save :: kpass=0
!  integer(ip)       :: kk,ii,mcomp
  !
  ! RESGS_TEM: Compute SGS residual
  !
  !if( IPARALL ) then
  !   nparr =  2
  !   parre => resgs_tem
  !   call par_operat(3_ip)
  !end if

  !
  ! Maximum subgrid scale residual
  !
  !
  ! Compute SGS residual
  !
  !if(kfl_sgsno_tem==1.or.kfl_sgsti_tem==1) then
  !   if(kfl_paral>=0) then
  !      nparr    =  2
  !      parre    => resgs_tem
  !      call par_operat(3_ip)
  !   end if
  !!   if(resgs_tem(2)>0.0_rp) resgs_tem(1)=sqrt(resgs_tem(1)/resgs_tem(2))      
  !end if
 !
  ! Write convergence
  !
  !if(kfl_paral<=0)&
  !     write(momod(modul)%lun_conve,'(110x,a,e15.7)') 'sgsres=',resgs_tem(1)


!!$  if( kfl_sgsco_tem /= 0 .and. IPARALL ) then
!!$     nparr =  1
!!$     parre => rmsgs_tem
!!$     call par_operat(2_ip)
!!$  end if
!!$  !
!!$  ! Subgrid scale statistics for convergence 
!!$  !
!!$  if( kfl_sgsco_tem /= 0 ) then
!!$
!!$     if( kpass==0 .and. kfl_rstar/=2 .and. INOTSLAVE ) then              
!!$        write(lun_stasg_tem,200) tosgs_tem,misgs_tem
!!$        write(lun_cvgsg_tem,300)
!!$     end if
!!$
!!$     if( IPARALL ) then
!!$        npari =  misgs_tem
!!$        parin => itsta_tem
!!$        call par_operat(3_ip)
!!$     end if
!!$
!!$     ii=0
!!$     do kk=1,misgs_tem
!!$        ii=ii+itsta_tem(kk)
!!$     end do
!!$     if(ii==0) ii=1
!!$
!!$     if( INOTSLAVE )&
!!$          write(lun_stasg_tem,201) &
!!$          (100.0_rp*real(itsta_tem(kk))/real(ii),kk=1,misgs_tem)
!!$
!!$     domcomp: do mcomp=misgs_tem,1,-1
!!$        if(itsta_tem(mcomp)/=0) exit domcomp
!!$     end do domcomp
!!$
!!$     if( IPARALL ) then
!!$        call pararr('SUM',0_ip,2*mcomp,resis_tem)
!!$     end if
!!$
!!$     if( INOTSLAVE ) then
!!$        do kk = 1,mcomp
!!$           if(resis_tem(2,kk)<1.0e-10) resis_tem(2,kk)=1.0_rp
!!$           resis_tem(1,kk)=sqrt(resis_tem(1,kk)/resis_tem(2,kk))
!!$           write(lun_cvgsg_tem,301) &
!!$                ittim,itcou,itinn(modul),kk,cutim,resis_tem(1,kk)
!!$        end do
!!$        flush(lun_stasg_tem)
!!$        flush(lun_cvgsg_tem)
!!$     end if
!!$
!!$  end if
!!$  kpass = 1 

  !----------------------------------------------------------------------
  !
  ! Formats
  !
  !----------------------------------------------------------------------

200 format(&
       & '# Percentage of the number of iterations to achieve the L2 tolerance= ',e12.6,&
       & ' with a maxiumum number of iterations= ',i3,/,&
       &      '#            1','             2','             3','             4',&
       &      '             5','             6','             7','             8',&
       &      '             9','            10','            11','            12',&
       &      '            13','            14','            15','           >16')
201 format(500(2x,e12.6))
300 format('# --| ALYA convergence '  ,/,&
       & '# --| Columns displayed:' ,/,&
       & '# --| 1. Time Step         2. Global Iteration  3. Inner Iteration   ',/,&
       & '# --| 4. Inner SGS         5. Current time      6. Velocity SGS      ',/,&
       & '# ','          1','          2','          3','          4',&
       &      '             5','             6') 
301 format(2x,4(2x,i9),3(2x,e12.6))

end subroutine tem_cvgsgs
