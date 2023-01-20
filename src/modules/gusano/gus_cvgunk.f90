!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Gusano
!> @{
!> @file    gus_cvgunk.f90
!> @author  houzeaux
!> @date    2020-10-22
!> @brief   Convergence
!> @details Compute and check convergence
!> @} 
!-----------------------------------------------------------------------

subroutine gus_cvgunk(itask)

  use def_parame
  use def_master
  use def_domain
  use def_gusano
  use mod_communications,       only : PAR_AVERAGE
  use mod_ker_detection,        only : ker_detection_doiter
  use mod_outfor,               only : outfor
  use mod_iofile,               only : iofile_flush_unit
  use mod_array_operations,     only : array_operations_residual_norm
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip), save       :: kpass=0
  real(rp),    save       :: cpuit_gus=0.0_rp,rigus=0.0_rp,prgus=0.0_rp
  real(rp),    save       :: rL2fl=0.0_rp,rL2pr=0.0_rp
  real(rp),    save       :: rL2pm=0.0_rp,rL2pc=0.0_rp
  real(rp)                :: time1
!  real(rp)                :: dtdtc
  character(21)           :: wnor1,wnor2,wnor3,wnor4
  real(rp)                :: vemin_gus,vemax_gus
  
  select case ( itask )

  case ( 0_ip )

     !-------------------------------------------------------------------
     !
     ! Timings and output
     !
     !-------------------------------------------------------------------

     wnor1 = 'L2 flow rate  '
     wnor2 = 'L2 pressure   '
     wnor3 = 'L2 velocity   '
     wnor4 = 'L2 pressure   '
     !
     ! Write convergence
     !
     if( INOTSLAVE ) then
        call cputim(time1)
        if( kpass == 0 .and. kfl_rstar /= 2 ) then
           write(momod(modul) % lun_conve,100) wnor1,wnor2
        end if
        if( kpass == 1 ) then
           time1 = time1 - cpuit_gus
        else
           time1 = time1 - cpu_initi
        end if
        !if( dtcri_gus /= 0.0_rp ) then
        !   dtdtc = dtime / dtcri_gus
        !else
        !   dtdtc = 0.0_rp
        !end if
        !dtdtc = dtcri_gus
        vemin_gus = 0.0_rp
        vemax_gus = 0.0_rp
        write(momod(modul) % lun_conve,101) &
             !   1     2            3     4     5     6     7     8
             ittim,itcou,itinn(modul),cutim,rL2fl,rL2pr,rL2pm,rL2pc,&
             !       9        10   
             vemin_gus,vemax_gus,&
             !  11    12    13
             rigus,prgus,time1

        call cputim(cpuit_gus)
        call iofile_flush_unit(momod(modul) % lun_conve)
     end if
     kpass=1     

  case ( ITASK_ENDINN )

     !-------------------------------------------------------------------
     !
     ! Check convergence of the inner iterations
     !
     !-------------------------------------------------------------------
     !
     ! RL2VE: Flow rate L2 residual
     ! RL2PR: Pressure L2 residual
     ! RL2SG: SGS L2 residual
     !
     if( solve(1) % kfl_block == 0 ) then
        rL2fl = array_operations_residual_norm(2_ip,2_ip,1_ip,unkno    ,flowr    ,0_ip  ,0_ip ,1_ip)
        rL2pr = array_operations_residual_norm(2_ip,2_ip,1_ip,unkno    ,press    ,1_ip  ,0_ip ,1_ip)
     else
        rL2fl = array_operations_residual_norm(2_ip,1_ip,1_ip,unkno    ,flowr    ,0_ip  ,0_ip ,1_ip)
        rL2pr = array_operations_residual_norm(2_ip,1_ip,1_ip,unkno    ,press    ,npoin ,0_ip ,1_ip)
     end if
     rL2pm = array_operations_residual_norm(2_ip,1_ip,1_ip,projm_gus,projm_gus,0_ip,npoin,1_ip)
     rL2pc = array_operations_residual_norm(2_ip,1_ip,1_ip,projc_gus,projc_gus,0_ip,npoin,1_ip)
     !
     ! KFL_GOITE_GUS: Check convergence
     !
     if( isnain(rL2fl)     .or.  isnain(rL2pr)     ) kfl_goite_gus = 0  ! NaN or +/- Inf
     if( rL2fl > 1.0e10_rp .or.  rL2pr > 1.0e10_rp ) kfl_goite_gus = 0  ! Huge
     if( rL2fl < cotol_gus .and. rL2pr < cotol_gus ) kfl_goite_gus = 0  ! Inner tolerance achieved
     if( itinn(modul) >= momod(modul) % miinn      ) kfl_goite_gus = 0  ! Number of iterations achieved
     
  case ( ITASK_ENDITE )
     !
     ! Check convergence of the outer iterations:
     ! || T(n,i,*) - T(n,i-1,*)|| / ||T(n,i,*)||
     !
     !resid_tem = array_operations_residual_norm(2_ip,1_ip,1_ip,flowr,flowr,0_ip,npoin,1_ip) 
     
  case( ITASK_ENDSTE )
     !
     ! Check residual of the time iterations:
     ! || T(n,*,*) - T(n-1,*,*)|| / ||T(n,*,*)||
     !
     rigus = array_operations_residual_norm(2_ip,1_ip,1_ip,flowr,flowr,0_ip,2_ip*npoin,1_ip)       
     prgus = array_operations_residual_norm(2_ip,1_ip,1_ip,press,press,0_ip,2_ip*npoin,1_ip)       
     if(rigus<=sstol_gus) then
        kfl_stead_gus = 1
        call outfor(28_ip,momod(modul)%lun_outpu,' ')
     end if
     if( kfl_timei_gus /= 0 ) kfl_stead_gus = 0
     
  end select
  
100 format('# --| ALYA convergence  ' ,/,&
       & '# --| Columns displayed:' ,/,&
       & '# --|  1. Time Step             2. Global Iteration       3. Inner Iteration   ',/,&
       & '# --|  4. Current time          5. ',a21,              '  6. ',a21,          ' ',/,& 
       & '# --|  7. Velocity SGS          8. SGS Max residual       9. Min. velocity     ',/,&
       & '# --| 10. Max. velocity        11. Min. tau              12. Max. tau          ',/,&
       & '# --| 13. Total mass           14. Relax. velocity       15. Relax. pressure   ',/,&
       & '# --| 16. ',a21,              '17. ',a21,              ' 18. Dt/Dt_critical    ',/,&
       & '# --| 19. Ass. ave cpu time    20. Ass. max cpu time     21. Sol. ave cpu time ',/,&
       & '# --| 22. Sol. max cpu time    23. Mom. Stead. res.      24. Cont. Stead. res. ',/,&
       & '# --| 25. Safety factor        26. Min. pressure         27. Max. pressure     ',/,&
       & '# --| 28. Linf Algeb. mom.     29. Linf Algeb. cont.     30. Elapsed CPU time  ',/,&
       & '# --| 31. SGS ave cpu time     32. SGS max cpu time      33. Bou. ave cpu time ',/,&
       & '# --| 34. Bou. max cpu time    35. Time residual         ',/,&
       & '# ','          1','          2','          3',&
       &      '             4','             5','             6','             7',&
       &      '             8','             9','            10','            11',&
       &      '            12','            13','            14','            15',&
       &      '            16','            17','            18','            19',&
       &      '            20','            21','            22','            23',&
       &      '            24','            25','            26','            27',&
       &      '            28','            29','            30','            31',&
       &      '            32','            33','            34') 
101 format(4x,i9,2x,i9,2x,i9,42(2x,e12.6))
  
end subroutine gus_cvgunk
