!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_reanut()
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_reanut
  ! NAME 
  !    tur_reanut
  ! DESCRIPTION
  !    This routine reads the numerical treatment 
  ! USED BY
  !    tur_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_solver
  use def_turbul
  use def_domain
  use mod_memchk
  use mod_ecoute, only :  ecoute
  implicit none
  integer(ip) :: ivari

  if( INOTSLAVE ) then
     !
     !  Initializations (defaults)
     !
     kfl_algor_tur = 1                                ! Algorithm: coupled(2)/uncoupled(1)
     kfl_clipp_tur = 0                                ! Clipping strategy            
     kfl_ellen_tur = 0                                ! Minimum element length
     kfl_normc_tur = 2                                ! L2 norm for convergence
     kfl_repro_tur = 0                                ! Res. projection not used 
     kfl_taust_tur = 1                                ! Tau calculation option
     kfl_shock_tur = 0                                ! Shock capturing off
     kfl_relax_tur = 1                                ! Constant relaxation
     kfl_tiacc_tur = 1                                ! First order time integ.
     kfl_tisch_tur = 1                                ! Time intgeration scheme
     kfl_walgo_tur = 1                                ! Wall distance algorithm: use equation
     kfl_weigh_tur = 1                                ! dT/dt is in the residual
     kfl_assem_tur = 0                                ! Assembly type 
     kfl_ortho_tur = 0                                ! ASGS
     kfl_limit_tur = 0                                ! No limiter
     kfl_produ_tur = 0                                ! Discontinuous production term
     kfl_meshi_tur = 0                                ! Mesh interpolator activation (OFF=0,ON=1)
     miinn_tur     = 1                                ! One internal iteration
     niter_tur     = 1                                ! One inner iteration
     neule_tur     = 0                                ! Number of Euler time steps
     kfl_sgsti_tur = 0                                ! Subscale time tracking
     kfl_sgsno_tur = 0                                ! Subscale non-linear tracking
     kfl_tibub_tur = 0                                ! Time integration of bubble
     kfl_sgsac_tur = 1                                ! SGS time accuracy

     staco_tur(1)  = 1.0_rp                           ! Diffusive term
     staco_tur(2)  = 1.0_rp                           ! Convective term
     staco_tur(3)  = 1.0_rp                           ! Reactive term
     shock_tur     = 0.0_rp                           ! SC parameter
     sstol_tur     = 1.0e-5_rp                        ! Steady-state tolerance
     safet_tur     = 1.0e10_rp                        ! Safety factor
     cotol_tur     = 1.0_rp                           ! Internal tolerance
     relax_tur     = 1.0_rp                           ! Relaxation factor
     safex_tur     = 1.0_rp                           ! Time function parameter for safety factor
     safma_tur     = 1.0e9_rp                         ! Maximum safety factor
     safeo_tur     = safet_tur                        ! Initial safety factor     
     saflo_tur     = safet_tur                        ! Global safety factor for local time step
     bemol_tur     = 0.0_rp                           ! Bemol     
     solve_sol     => solve(1:)
     kfl_lmaxi_tur =  0

     clipfac_tur= 1.0e-6_rp

     ivari         = 1                                ! Variable used toread solver's parameter
     !
     ! Reach the section
     !
     call ecoute('tur_reanut')
     do while(words(1)/='NUMER')
        call ecoute('tur_reanut')
     end do
     !
     ! Begin to read data
     !
     do while(words(1)/='ENDNU')
        call ecoute('tur_reanut')

        if(words(1)=='TAUST') then

           call reatau(kfl_taust_tur)

        else if(words(1)=='TRACK') then

           if(exists('TIME ')) kfl_sgsti_tur = 1
           if(exists('NONLI')) then
              kfl_sgsno_tur = 1
           end if
           kfl_sgsac_tur = getint('ORDER',1_ip,'#Time integration order')
           if(exists('FIRST')) kfl_sgsac_tur = 1 
           if(exists('SECON')) kfl_sgsac_tur = 2

        else if( words(1) == 'STRAT' .or. words(1) == 'STABI' ) then
           if( words(2) == 'OFF  '.or. words(2) == 'GALER' ) then
              kfl_ortho_tur = -2
           else if( words(2) == 'SU   ' .or. words(2) == 'FIRST' ) then
              kfl_ortho_tur = -1
           else if( words(2) == 'ASGS ' ) then
              kfl_ortho_tur =  0
           else if( words(2) == 'FULLO' ) then
              kfl_ortho_tur =  1
           else if( words(2) == 'OSS  ' .or. words(2) == 'AOSS ' ) then
              kfl_ortho_tur =  2
              if( exists('NOLIM') ) kfl_limit_tur =  0  ! No limiter
              if( exists('SOTO ') ) kfl_limit_tur =  1  ! Soto limiter
              if( exists('DIFFU') ) kfl_limit_tur =  2  ! Very diffusive limiter
              if( exists('FIRST') ) kfl_limit_tur = -1  ! First order
           else if( words(2) == 'AROSS' ) then
              kfl_ortho_tur =  3
              if( exists('NOLIM') ) kfl_limit_tur =  0  ! No limiter
              if( exists('SOTO ') ) kfl_limit_tur =  1  ! Soto limiter
              if( exists('DIFFU') ) kfl_limit_tur =  2  ! Very diffusive limiter
              if( exists('FIRST') ) kfl_limit_tur = -1  ! First order
           else if( words(2) == 'SUPG ' ) then
              kfl_ortho_tur =  -3
           else if( words(2) == 'BUBBL' ) then
              kfl_ortho_tur =  -4
              if( words(3) == 'TIMET' ) kfl_tibub_tur = 1
           end if

        else if(words(1)=='BEMOL') then
           bemol_tur=getrea('BEMOL',0.0_rp,'#Bemol of convective term')

        else if(words(1)=='PRODU') then
           if( words(2) == 'SMOOT') then
              kfl_produ_tur = 1
           end if

        else if(words(1)=='ELEME') then
           call realen(kfl_ellen_tur)

        else if(words(1)=='TYPEO') then

           if(words(2)=='RESID') kfl_repro_tur = 1

        else if(words(1)=='SHOCK') then
           if(exists('ISOTR').or.exists('ON   ')) then
              kfl_shock_tur = 1
              if(exists('VALUE')) &
                   shock_tur     = getrea('VALUE',0.0_rp,'#Shock capturing parameter')
           else if(exists('ANISO')) then
              kfl_shock_tur = 2
              if(exists('VALUE')) &
                   shock_tur     = getrea('VALUE',0.0_rp,'#Shock capturing parameter')
           end if
        else if(words(1)=='TEMPO') then
           if(exists('GALER')) kfl_weigh_tur = 0
           if(exists('ALL  ')) kfl_weigh_tur = 1

        else if(words(1)=='ASSEM') then
           if(words(2)=='CELL ') then
              kfl_assem_tur=2
           end if
           
        else if(words(1)=='TIMEA') then
           kfl_tiacc_tur = int(param(1))
           if(kfl_timei_tur==0) kfl_tiacc_tur = 1
           neule_tur= getint('EULER',0_ip,'#EULER TIME STEPS')

        else if(words(1)=='TIMEI') then
           if(exists('TRAPE')) kfl_tisch_tur=1
           if(exists('BDF  ')) kfl_tisch_tur=2
           kfl_tiacc_tur = getint('ORDER',1_ip,'#Time integration order')
           neule_tur     = getint('EULER',0_ip,'#EULER TIME STEPS')

        else if(words(1)=='SAFET') then
           safet_tur = param(1)
           safeo_tur=  safet_tur ! initial safety factor
           
           safex_tur  = getrea('EXPON',1.0_rp,'#Safety Factor Exponential time function')
           
           safma_tur  = getrea('MAXIM',1.0e9_rp,'#Maximum safety factor function')

           saflo_tur  = getrea('MINGL',safet_tur,'#Minimum global safety factor when using local time step. Fixes minimum time step over elements')    

        else if(words(1)=='STEAD') then
           sstol_tur = param(1)

        else if(words(1)=='NORMO') then
           if(exists('L1   ')) then
              kfl_normc_tur = 1
           else if(exists('L2   ')) then
              kfl_normc_tur = 2
           else if(exists('L-inf')) then
              kfl_normc_tur = 0
           else if(exists('ALGEB')) then
              kfl_normc_tur = 3
           end if

        else if(words(1)=='MAXIM') then
           miinn_tur = int(param(1))
           
        else if(words(1)=='INNER') then
           niter_tur = int(param(1))

        else if(words(1)=='CONVE') then
           cotol_tur = param(1)

        else if(words(1)=='RELAX') then
           relax_tur = param(1)
           if(exists('CONST')) then
              kfl_relax_tur=1
           else if(exists('AITKE')) then
              kfl_relax_tur=2
           end if

        else if(words(1)=='ALGOR') then
           if(words(2)=='COUPL') then
              kfl_algor_tur=min(2_ip,nturb_tur)
           else
              kfl_algor_tur=1
           end if

        else if( words(1) == 'MESHI' ) then
           if(exists('ON   ') ) kfl_meshi_tur = 1

        else if( words(1) == 'TURBU' ) then
           ivari = 1

        else if( words(1) == 'LMAXI' ) then
           ivari = 3
           kfl_lmaxi_tur =  1

        else if(words(1)=='ALGEB') then
           solve_sol => solve(ivari:)
           call reasol(1_ip)
           solve_sol => solve(1:)

        else if(words(1)=='PRECO') then 
           solve_sol => solve(ivari:)
           call reasol(2_ip)
           solve_sol => solve(1:)

        else if(words(1)=='CLIPP') then
           if(words(2)=='NONE ') then
              kfl_clipp_tur=-2
           else if(words(2)=='OFF  ') then
              kfl_clipp_tur=-1
           else if(words(2)=='ALMOS') then
              kfl_clipp_tur=0
           else if(words(2)=='ZERO ') then
              kfl_clipp_tur=1
           else if(words(2)=='CHARA') then
              kfl_clipp_tur=2
           else if(words(2)=='FIXED') then
              kfl_clipp_tur=3
           else if(words(2)=='ABSOL') then
              kfl_clipp_tur=4
           else if(words(2)=='LASTV') then
              kfl_clipp_tur=5
           else if(words(2)=='PRESC') then
              kfl_clipp_tur=6
              clipfac_tur  = getrea('FACTO',1.0e-6_rp,'#Clipping factor percentage')
           end if

        end if
     end do
     !
     ! Copy solver data
     !
     if( kfl_algor_tur == 1 ) then
        solve_sol => solve(1:)
        if( nturb_tur >= 2 ) call solcpy(1_ip,2_ip)
        if( nturb_tur >= 3 ) call solcpy(1_ip,3_ip)
        if( nturb_tur >= 4 ) call solcpy(1_ip,4_ip)
     end if

  end if

end subroutine tur_reanut
