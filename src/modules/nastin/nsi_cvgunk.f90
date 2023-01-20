!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_cvgunk(itask)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_cvgunk
  ! NAME 
  !    nsi_cvgunk
  ! DESCRIPTION
  !    This routine performs several convergence checks for NASTIN
  ! USES
  !    nsi_endite (itask=1,2)
  !    nsi_endste (itask=3)
  ! USED BY
  !    Nastin
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_nastin
  use mod_nsi_schur_complement, only : nsi_momentum_continuity_residuals
  use mod_communications,       only : PAR_AVERAGE
  use mod_ker_detection,        only : ker_detection_doiter
  use mod_outfor,               only : outfor
  use mod_iofile,               only : iofile_flush_unit
  use mod_array_operations,     only : array_operations_residual_norm
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip), save       :: kpass=0,jpass=0
  real(rp),    save       :: cpuit_nsi=0.0_rp,rinsi,prnsi,resti=1.0_rp
  real(rp),    save       :: xmass
  real(rp),    save       :: rLive=0.0_rp,rLipr=0.0_rp
  real(rp)                :: time1,dtdtc
  character(21)           :: wnor1,wnor2,wnor3,wnor4
  real(rp)                :: cpu_ass_max
  real(rp)                :: cpu_ass_ave
  real(rp)                :: cpu_bou_max
  real(rp)                :: cpu_bou_ave
  real(rp)                :: cpu_sol_max
  real(rp)                :: cpu_sol_ave
  real(rp)                :: cpu_sgs_max
  real(rp)                :: cpu_sgs_ave

  select case ( itask )

  case ( 0_ip )

     !-------------------------------------------------------------------
     !
     ! Timings and output
     !
     !-------------------------------------------------------------------

     wnor1 = 'Linf velocity '
     wnor2 = 'Linf pressure '
     wnor3 = 'L2 Algeb. mom.'
     wnor4 = 'L2 Algeb. con.'
     !
     ! Min and max: veloc, press, tau
     !
     call nsi_minmax(&
          cpu_ass_ave,cpu_ass_max,cpu_bou_ave,cpu_bou_max,&
          cpu_sol_ave,cpu_sol_max,cpu_sgs_ave,cpu_sgs_max)
     !
     ! Write convergence
     !
     if( INOTSLAVE ) then
        call cputim(time1)
        if( kpass == 0 .and. kfl_rstar /= 2 ) then
           write(momod(modul) % lun_conve,100) wnor1,wnor2,wnor3,wnor4
        end if
        if( kpass == 1 ) then
           time1 = time1 - cpuit_nsi
        else
           time1 = time1 - cpu_initi
        end if
        if( dtcri_nsi /= 0.0_rp ) then
           dtdtc = dtime / dtcri_nsi
        else
           dtdtc = 0.0_rp
        end if
        dtdtc = dtcri_nsi

        write(momod(modul) % lun_conve,101) &
             !   1     2            3     4     5     6            7
             ittim,itcou,itinn(modul),cutim,rLive,rLipr,resgs_nsi(1),&
             !       8         9        10        11        12    13        14   
             rmsgs_nsi,vemin_nsi,vemax_nsi,tamin_nsi,tamax_nsi,xmass,relax_nsi,&
             !      15    16    17    18 
             relap_nsi,rinsi,prnsi,dtdtc,&
             !        19          20          21          22
             cpu_ass_ave,cpu_ass_max,cpu_sol_ave,cpu_sol_max,&
             !         23           24        25        26        27 
             resss_nsi(1),resss_nsi(2),safet_nsi,prmin_nsi,prmax_nsi,&
             !         28           29    30          31          32
             reinf_nsi(1),reinf_nsi(2),time1,cpu_sgs_ave,cpu_sgs_max,&
             !        34          34    35
             cpu_bou_ave,cpu_bou_max,resti

        call cputim(cpuit_nsi)
        call iofile_flush_unit(momod(modul) % lun_conve)
     end if
     kpass=1     

  case ( ITASK_ENDINN )

     !-------------------------------------------------------------------
     !
     ! Check convergence of the inner iterations
     ! RLIVE: Velocity Linf residual
     ! RLIPR: Pressure Linf residual
     ! RINSI: Velocity algebraic residual
     ! PRNSI: Pressure algebraic residual
     !
     !-------------------------------------------------------------------

     if( NSI_MONOLITHIC ) then
        rLive = array_operations_residual_norm(0_ip,ndime+1,ndime,unkno,veloc    ,0_ip     ,0_ip,ndime)      
        rLipr = array_operations_residual_norm(0_ip,ndime+1, 1_ip,unkno,unk2n_nsi,ndime    ,0_ip,1_ip)      
     else if( NSI_FRACTIONAL_STEP .or. NSI_SEMI_IMPLICIT ) then
        rLive = array_operations_residual_norm(0_ip,ndime  ,ndime,unkno,veloc    ,0_ip     ,ndime*npoin,ndime)
        rLipr = array_operations_residual_norm(0_ip, 1_ip  , 1_ip,unkno,unk2n_nsi,ndbgs_nsi,npoin,1_ip)
     else
        rLive = array_operations_residual_norm(0_ip,ndime  ,ndime,unkno,veloc    ,0_ip     ,0_ip,ndime)   
        rLipr = array_operations_residual_norm(0_ip,1_ip   , 1_ip,unkno,unk2n_nsi,ndbgs_nsi,0_ip,1_ip)   
     end if
     
     if( kfl_normc_nsi == 4 ) then                                                                    ! Error w/r reference solution
        call nsi_refere()                                                                             ! Reference solution
        rinsi = difve_nsi
        prnsi = difpr_nsi
     else     
        if( NSI_MONOLITHIC ) then
           call nsi_residual_monolithic(amatr,rhsid,veloc,press,resin_nsi(1),resin_nsi(2))
        else if( NSI_SCHUR_COMPLEMENT ) then
           call nsi_momentum_continuity_residuals(1_ip)               
        else if( NSI_FRACTIONAL_STEP .or. NSI_SEMI_IMPLICIT ) then
           resin_nsi(1) = array_operations_residual_norm(2_ip,ndime  ,ndime,unkno,veloc    ,0_ip     ,ndime*npoin,ndime)
           resin_nsi(2) = array_operations_residual_norm(2_ip, 1_ip  , 1_ip,unkno,unk2n_nsi,ndbgs_nsi,npoin,1_ip)
        end if
        rinsi = resin_nsi(1)
        prnsi = resin_nsi(2)
     end if
     !
     ! KFL_GOITE_NSI: Check convergence
     !
     if( isnain(rinsi)  .or. isnain(prnsi)  )  kfl_goite_nsi = 0                               ! NaN or +/- Inf
     if( rinsi > 1.0e10_rp .or. prnsi > 1.0e10_rp )  kfl_goite_nsi = 0                         ! Huge
     if( (rinsi < cotol_nsi .and. prnsi < cotol_nsi ) .or. &                                   ! Inner tolerance achieved
          itinn(modul) >= momod(modul) % miinn ) kfl_goite_nsi = 0                             ! Number of iterations achieved
     if( kfl_refer_nsi /= 0 ) then
        if( isnain(difve_nsi)    .or. isnain(difpr_nsi)    ) kfl_goite_nsi = 0                 ! NaN or +/- Inf        
        if( difve_nsi > 1.0e7_rp .or. difpr_nsi > 1.0e7_rp ) kfl_goite_nsi = 0                 ! Huge    
     end if
     !
     ! Detection of non-convergence
     !
     call ker_detection_doiter(ndime,ndime,unkno,veloc,one,one,ndime,'INNER_NOT_CONVERGED','VELOC','NO_COMMENT')
     !
     ! XMASS: Compute mass and normalize it
     !
     xmass = 0.0_rp
     !call nsi_dommas(xmass)

  case ( ITASK_ENDITE )

     !-------------------------------------------------------------------
     !
     ! Check convergence of the outer iterations:
     ! RESID_NSI: velocity residual
     ! RESIP_NSI: pressure residual 
     !
     !-------------------------------------------------------------------

     if( kfl_normc_nsi == 4 ) then                                       ! L2 Error w/r reference solution
        resid_nsi = difve_nsi
        resip_nsi = difpr_nsi
     else                                                                ! L2 norm: ||u^{n^+1,i+1}-u^{n+1,i-1}||/||u^{n+1,i+1}||
        resid_nsi = array_operations_residual_norm(2_ip,ndime,ndime,veloc,    veloc    ,0_ip,ndime*npoin,ndime)
        resip_nsi = array_operations_residual_norm(2_ip, 1_ip, 1_ip,unk2n_nsi,unk2n_nsi,0_ip,      npoin, 1_ip)
     end if

  case ( ITASK_ENDSTE ) 

     !-------------------------------------------------------------------
     !
     ! Check residual of the time iterations
     ! RESTI: time residual measured only using velocity
     !
     !-------------------------------------------------------------------

     if( kfl_normc_nsi == 4 .and. kfl_timei_nsi /= 0 ) then              ! L2 wrt Reference solution
        resti = max(difve_nsi,difpr_nsi)
     else if( kfl_timei_nsi /= 0 ) then                                  ! L2 residual: ||u^{n+1}-u^n||/||u^{n+1}||
        resti = array_operations_residual_norm(2_ip,ndime,ndime,veloc,veloc,0_ip,2_ip*ndime*npoin,ndime)
     end if
     !
     ! Some check to stop the code if things are going bad
     !
     if( isnain(resti)      ) kfl_stead_nsi = 1                          ! NaN or +/- Inf
     if( resti > 1.0e10_rp  ) kfl_stead_nsi = 1                          ! Huge
     if( kfl_refer_nsi /= 0 ) then
        if( isnain(difve_nsi) .or. isnain(difpr_nsi) ) kfl_stead_nsi = 1 ! NaN or +/- Inf        
     end if
     !
     ! Perform at least two time step
     ! Example: coupling with temperature where initial temperature is constant
     !          and flow confined: it starts to move at second iteration.
     !
     if( resti <= sstol_nsi .and. ittim >= 2 ) then
        kfl_stead_nsi = 1
        call outfor(28_ip,momod(modul) % lun_outpu,' ')
     end if
     !
     ! Low-Mach model
     !
     if( INOTSLAVE .and. kfl_regim_nsi==3 ) then
        if( jpass == 0 ) then
           jpass = 1
           write(lun_lmach_nsi,400)
        end if
        write(lun_lmach_nsi,401) cutim,prthe(1),prthe(1)/prthe(4),dpthe,xmass_nsi
     end if

  end select

  !----------------------------------------------------------------------
  !
  ! Formats
  !
  !----------------------------------------------------------------------

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
400 format('# --| ALYA Low-Mach model variables '  ,/,&
       & '# --| Columns displayed:' ,/,&
       & '# --| 1. Time Step         2. Therm. pres. p0   3. p0/p0^0            ',/,&
       & '# --| 4. dp/dt             5. Total mass        ',/,&
       & '# ','             1','             2','             3','             4',&
       &      '             5') 
401 format(2x,10(2x,es16.8e3))

end subroutine nsi_cvgunk
