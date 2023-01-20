!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_reanut
  !-----------------------------------------------------------------------
  !****f* Levels/lev_reanut 
  ! NAME 
  !    lev_reanut
  ! DESCRIPTION
  !    This routine reads the numerical treatment for LEVELS module
  ! USES
  !    ecoute
  ! USED BY
  !    lev_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_solver
  use def_levels
  use def_domain
  use mod_memchk
  use mod_ecoute, only :  ecoute
  implicit none
  integer :: iprob,iredi

  if( INOTSLAVE ) then
     ! 
     !  Initializations (defaults)
     !
     kfl_timet_lev = 1                                ! Explicit(1)/Implicit(2)

     kfl_tisch_lev = 2                                ! BDF
     kfl_timco_lev = 1                                ! Adapt time step just for the LS
     kfl_tiacc_lev = 1                                ! First order time integ.
     kfl_normc_lev = 0  
     kfl_ellen_lev = 0                                ! Minimum element length
     kfl_zonal_lev = 0                                ! Redistanzation by zones (default all domain )
     kfl_geodi_lev = 1                                !Cristobal's geometrical distance method by default
     neule_lev     = 0                                ! Number of Euler time steps  
     miinn_lev     = 1
     inred_lev     = 0                                ! Initial Redistanciation
     nfred_lev     = 10                               ! Redistanciation frequency
     tyred_lev     = 0                                ! Redistanciation type
     nstre_lev     = 0                                ! If /=0 indicates the step to stop redistancing, it is used to debug error with redistance
     kfl_locre_lev = 0                                ! Local redistanciation, in a small layer (thicl*10). In the rest step over with the value prior to redistance
     nbitr_lev     = 1                                ! Redistanciation  equation iteration number
     kfl_corvo_lev = 0                                ! Volume correction throug
     safet_lev     = 1.0_rp                           ! Safety factor
     sstol_lev     = 1.0e-8_rp                        ! Steady state tolerance
     cotol_lev     = 0.1_rp                           ! Internal tolerance
     cpuit_lev     = 0.0_rp                           ! CPU time per iteration
     supgp_lev     = 1.0_rp                           ! SUPG stabililization 
     solve_sol     => solve                           ! Solver type
     !
     ! Reach the section
     !
     iprob = 1
     iredi = 0
     call ecoute('lev_reanut')
     do while(words(1)/='NUMER')
        call ecoute('lev_reanut')
     end do
     !
     ! Begin to read data
     !
     do while(words(1)/='ENDNU')
        call ecoute('lev_reanut')

        if(words(1)=='TIMET') then
           if(words(2)=='EXPLI') then
              kfl_timet_lev=1
           else if(words(2)=='IMPLI') then
              kfl_timet_lev=2
           end if

        else if(words(1)=='ELEME') then
           call realen(kfl_ellen_lev)

        else if(words(1)=='TIMEI') then
           if(exists('TRAPE')) kfl_tisch_lev=1
           kfl_tiacc_lev = getint('ORDER',1_ip,'#Time integration order')
           if(exists('RUNGE')) kfl_tisch_lev=4
           kfl_tiacc_lev = getint('ORDER',1_ip,'#Time integration order')
           if(exists('BDF  ')) kfl_tisch_lev=2
           kfl_tiacc_lev = getint('ORDER',1_ip,'#Time integration order')
           neule_lev     = getint('EULER',0_ip,'#EULER TIME STEPS')
           kfl_timco_lev = getint('TIADA',0_ip,'#Adapt time step')

        else if(words(1)=='TIMEA') then
           kfl_tiacc_lev = int(param(1))
           neule_lev = getint('EULER',0_ip,'#EULER TIME STEPS')

        else if(words(1)=='SAFET') then
           safet_lev = param(1)

        else if(words(1)=='STEAD') then
           sstol_lev = param(1)

        else if(words(1)=='MAXIM') then
           miinn_lev = int(param(1))

        else if(words(1)=='CONVE') then
           cotol_lev = param(1)

        else if(words(1)=='NORMO') then
           if(exists('L2   ')) then
              kfl_normc_lev = 2
           else if(exists('L1   ')) then
              kfl_normc_lev = 1
           else if(exists('L-inf')) then
              kfl_normc_lev = 0
           end if

        else if(words(1)=='GENER') then
           iprob = 5
           
        else if(words(1)=='EQUAT') then
           if(     words(2)=='REDIS') then
              iprob = 2
           else if(words(2)=='MESHD') then
              iprob = 4
           else if(words(2)=='GENER') then
              iprob = 5
           else 
              iprob = 1
           end if
           solve_sol => solve(iprob:)

        else if(words(1)=='ALGEB') then
           solve_sol => solve(iprob:)
           call reasol(1_ip)
           iprob = 1

        else if(words(1)=='PRECO') then 
           solve_sol => solve(iprob:)
           call reasol(2_ip)
           iprob = 1
           
        else if(words(1)=='POISS') then
           solve_sol => solve(3:)     
           call ecoute('lev_reanut')
           do while(words(1)/='ENDPO')
              if(words(1)=='ALGEB') then
                 call reasol(1_ip)
              else if(words(1)=='PRECO') then 
                 call reasol(2_ip)
              end if
              call ecoute('lev_reanut')
           end do
           solve_sol => solve(1:)

        else if(words(1)=='SUSSM') then
           iprob = 2

        else if(words(1)=='REDIS') then
           inred_lev     = getint('INITI', 0_ip,'# INITIAL REDISTANCIATION ')
           nfred_lev     = getint('FREQU',10_ip,'#REDISTANCIATION FREQUENCY')
           tyred_lev     = getint('TYPER', 1_ip,'#REDISTANCIATION TYPE')
           nstre_lev     = getint('STOPA', 0_ip,'#STOP REDISTANCIATION AT STEP')
           nbitr_lev     = getint('NBITE', 1_ip,'#REDISTANCIATION EQUATION ITERATION NUMBER')
           kfl_zonal_lev = getint('ZONAL', 0_ip,'#REDISTANCIATION BY ZONES')
           if( exists('SAMES' ) ) iredi = 1
           if( exists('LOCAL' ) ) kfl_locre_lev = 1

        else if(words(1)=='GEOME') then !Type of geometrical distance
           if(words(2)  /='CRIST') kfl_geodi_lev = 2
           
        else if(words(1)=='VOLUM') then
           if(words(2)=='ON   ') kfl_corvo_lev=1
           if(words(2)=='SMOLI') kfl_corvo_lev=1
           if(words(2)=='LOEHN') kfl_corvo_lev=2

        else if(words(1)=='SOLVE') then
           solve(1)%miter= int(param(1))

        else if(words(1)=='TOLER') then
           solve(1)%solco = param(1)

        else if(words(1)=='KRYLO') then
           solve(1)%nkryd = int(param(1))
 
        end if

     end do
     !
     ! Copy level set solver to redistantiation solver if it has not been defined
     !
     if( iredi == 1 .and. solve(2)%kfl_algso==-999 ) then
        solve_sol => solve(1:)
        call solcpy(1_ip,2_ip)
     end if
     !
     ! temporary solution to avoid SOLDEF: NO SOLVER HAS BEEN DEFINED !!! search better solution
     !    
     if(solve(3)%kfl_algso==-999) then  
        solve_sol => solve(1:)
        call solcpy(1_ip,3_ip)
     end if

  end if

end subroutine lev_reanut
