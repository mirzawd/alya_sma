!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_reanut
  !-----------------------------------------------------------------------
  !****f* Temper/tem_reanut
  ! NAME 
  !    tem_reanut
  ! DESCRIPTION
  !    This routine reads the numerical treatment for TEMPER module
  ! USES
  !    ecoute
  ! USED BY
  !    tem_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_solver
  use def_temper
  use def_domain
  use mod_ecoute, only :  ecoute
  implicit none
  integer(ip) :: ivari

  if( INOTSLAVE ) then
     ! 
     !  Initializations (defaults)
     !
     kfl_dttyp_tem = 1                                ! Module time step strategy
     kfl_ellen_tem = 0                                ! Minimum element length
     kfl_sgsti_tem = 0                                ! Subgrid scale time tracking
     kfl_sgsno_tem = 0                                ! Subgrid scale non-lniear tracking
     kfl_taust_tem = 1                                ! Tau calculation option
     kfl_ortho_tem = 0                                ! Orthogonal SGS
     kfl_limit_tem = 0                                ! Lmiter
     kfl_shock_tem = 0                                ! Shock capturing off
     kfl_tiacc_tem = 1                                ! First order time integ.
     kfl_tibub_tem = 0                                ! Time integration of Bubble function
     kfl_assem_tem = 0                                ! Assembly type
     kfl_negat_tem = 0                                ! Strategy for negative overshoots
     kfl_posit_tem = 0                                ! Strategy for positive overshoots
     kfl_plepp_tem = -1                               ! PLE activation: =-1 OFF, =0 TEMPER, =3 LOWMACH, =4 ENTHALPY
     neule_tem     = 0                                ! Number of Euler time steps
     kfl_tisch_tem = 1                                ! Trapezoidal rule
     kfl_normc_tem = 2                                ! L2 norm for convergence
     miinn_tem     = 1                                ! One internal iteration
     misgs_tem     = 1                                ! Max # of SGS iterations
     kfl_sgsac_tem = 1
     kfl_sgsli_tem = 1                                ! SGS convection linearization PICARD
     kfl_meshi_tem = 0                                ! Mesh interpolator activation (OFF=0,ON=1)
     kfl_discr_tem = 0                                ! Discretization method (FE=1,FV=1)
     kfl_explicit_tem = 0                             ! explicit time integration
     kfl_rhs_scal_tem = 0                             ! Scaling PDE by density
     kfl_entropy_tem = 0                              ! entropy viscosity stablization method
     kfl_disable_entpre = 0                           ! Disables temporal term of entropy function along with RK prediction
     !
     staco_tem(1)  = 1.0_rp                           ! Diffusive term
     staco_tem(2)  = 1.0_rp                           ! Convective term
     staco_tem(3)  = 1.0_rp                           ! Reactive term
     shock_tem     = 0.0_rp                           ! SC parameter
     safet_tem     = 1.0e10_rp                        ! Safety factor
     source_safet_tem = 1.0_rp                        ! Source term safety factor
     sstol_tem     = 1.0e-5_rp                        ! Steady-state tolerance
     cotol_tem     = 0.1_rp                           ! Internal tolerance
     relax_tem     = 1.0_rp                           ! Relaxation factor
     bemol_tem     = 0.0_rp                           ! Bemol (convectice term)
     relsg_tem     = 1.0_rp                           ! Relaxation parameter of subgrid scale
     tosgs_tem     = 0.01_rp                          ! Subgrid scale tolerance
     solve_sol     => solve                           ! Solver type

     ivari        =  1                                ! Current solver
     !
     ! Reach the section
     !
     call ecoute('tem_reanut')
     do while(words(1)/='NUMER')
        call ecoute('tem_reanut')
     end do
     !
     ! Begin to read data
     !
     do while(words(1)/='ENDNU')
        call ecoute('tem_reanut')

        if(words(1)=='TRACK') then
           if(exists('TIME ')) kfl_sgsti_tem=1
           if(exists('NONLI')) then
              kfl_sgsno_tem=1
              misgs_tem=getint('ITERA',1_ip,   '#Subgrid scale iterations')
              tosgs_tem=getrea('TOLER',0.01_rp,'#Subgrid scale Tolerance')
              relsg_tem=getrea('RELAX',1.0_rp, '#Subgrid scale Relaxation') 
           end if
           kfl_sgsac_tem = getint('ORDER',1_ip,'#Time integration order')
           if(exists('FIRST')) kfl_sgsac_tem = 1 
           if(exists('SECON')) kfl_sgsac_tem = 2
           if(exists('PICAR')) kfl_sgsli_tem = 1
           if(exists('NEWTO')) kfl_sgsli_tem = 2

        else if(words(1)=='TAUST') then
           call reatau(kfl_taust_tem)       
           
        else if( words(1) == 'DISCR' ) then
           if( words(2) == 'FE   ' ) then
              kfl_discr_tem = 0
           else if( words(2) == 'FV   ' ) then
              kfl_discr_tem = 1
           end if

        else if(words(1)=='STRAT' .or. words(1) == 'STABI' ) then

           if( words(2) == 'SU   ' .or. words(2) == 'FIRST' ) then
              kfl_ortho_tem = -1
           else if( words(2) == 'ASGS ' ) then
              kfl_ortho_tem =  0
           else if( words(2) == 'FULLO' ) then
              kfl_ortho_tem =  1
           else if( words(2) == 'OSS  ' .or. words(2) == 'AOSS ' ) then
              kfl_ortho_tem =  2
              if( exists('NOLIM') ) kfl_limit_tem =  0  ! No limiter
              if( exists('SOTO ') ) kfl_limit_tem =  1  ! Soto limiter
              if( exists('DIFFU') ) kfl_limit_tem =  2  ! Very diffusive limiter
              if( exists('FIRST') ) kfl_limit_tem = -1  ! First order
           else if( words(2) == 'AROSS' ) then
              kfl_ortho_tem =  3
              if( exists('NOLIM') ) kfl_limit_tem =  0  ! No limiter
              if( exists('SOTO ') ) kfl_limit_tem =  1  ! Soto limiter
              if( exists('DIFFU') ) kfl_limit_tem =  2  ! Very diffusive limiter
              if( exists('FIRST') ) kfl_limit_tem = -1  ! First order
           else if ( words(2) == 'ENTRO') then
              if (words(3) == 'ENTPR') then
                 if (words(4) == 'OFF') then
                    kfl_disable_entpre = 1
                 else if (words(4) == 'ON') then
                    kfl_disable_entpre = 0
                 end if
              end if
              kfl_ortho_tem = -2
              kfl_entropy_tem = 1
           else if( words(2) == 'OFF  ' .or.  words(2) == 'GALER' ) then
              kfl_ortho_tem =  -2
           else if( words(2) == 'MARGA' ) then
              kfl_ortho_tem =  4
           else if( words(2) == 'SUPG ' ) then
              kfl_ortho_tem =  -3
           else if( words(2) == 'BUBBL' ) then
              kfl_ortho_tem =  -4
              if( words(3) == 'TIMET' ) kfl_tibub_tem = 1
           else if(words(2)=='CONST' ) then ! reads stab constants
              staco_tem(1) = param(2)
              staco_tem(2) = param(3)
              staco_tem(3) = param(4)              
           end if

        else if(words(1)=='TIMES') then
           if(words(2)=='LOCAL') then
              kfl_dttyp_tem=1
           end if

        else if(words(1)=='ASSEM') then
           if(words(2)=='ADR  ') then
              kfl_assem_tem=1
           else if(words(2)=='CELL ') then
              kfl_assem_tem=2
           end if

        else if(words(1)=='ELEME') then
           call realen(kfl_ellen_tem)
           if(words(2)=='NEW  ') kfl_ellen_tem=-1

        else if(words(1)=='SHOCK') then
           if(exists('ISOTR').or.exists('ON   ')) then
              kfl_shock_tem = 1
              shock_tem     = getrea('VALUE',0.0_rp,'#Shock capturing parameter')
           else if(exists('ANISO')) then
              kfl_shock_tem = 2
              shock_tem     = getrea('VALUE',0.0_rp,'#Shock capturing parameter')
           end if

        else if(words(1)=='TIMEI') then
           if(exists('TRAPE')) kfl_tisch_tem=1
           if(exists('BDF  ')) kfl_tisch_tem=2
           if(exists('ADAMS')) then
              kfl_tisch_tem = 3
              kfl_explicit_tem = 1
           end if
           if(exists('RUNGE')) then
              kfl_tisch_tem = 4
              kfl_explicit_tem = 1
           end if
           kfl_tiacc_tem = getint('ORDER',1_ip,'#Time integration order')
           neule_tem     = getint('EULER',0_ip,'#EULER TIME STEPS')

        else if(words(1)=='TIMEA') then
           kfl_tiacc_tem = int(param(1))
           if(kfl_timei_tem==0) kfl_tiacc_tem = 1
           neule_tem = getint('EULER',0_ip,'#EULER TIME STEPS')

        else if(words(1)=='SAFET') then
           safet_tem = param(1)
           if(exists('SOURC')) source_safet_tem = getrea('SOURC',1.0_rp,'#Source term safety factor')

        else if(words(1)=='STEAD') then
           sstol_tem = param(1)

        else if(words(1)=='RHSSC') then
           if(words(2)=='ON   ') kfl_rhs_scal_tem = 1

        else if(words(1)=='NORMO') then
           if(exists('L1   ')) then
              kfl_normc_tem = 1
           else if(exists('L2   ')) then
              kfl_normc_tem = 2
           else if(exists('L-inf')) then
              kfl_normc_tem = 0
           else if(exists('ALGEB')) then
              kfl_normc_tem = 3
           end if

        else if(words(1)=='MAXIM') then
           miinn_tem = int(param(1))

        else if(words(1)=='CONVE') then
           cotol_tem = param(1)

        else if(words(1)=='RELAX') then
           relax_tem = param(1)

        else if(words(1)=='BEMOL') then
           bemol_tem=getrea('BEMOL',0.0_rp,'#Bemol of convective term')

        else if( words(1) == 'MESHI' ) then
           if(exists('ON   ') ) kfl_meshi_tem = 1

        else if(words(1)=='CONSI') then
           ivari = 2
           
        else if(words(1)=='TEMPE') then
           ivari = 1
           
        else if(words(1)=='ALGEB') then
           solve_sol => solve(ivari:)
           call reasol(1_ip)

        else if(words(1)=='PRECO') then
           solve_sol => solve(ivari:)
           call reasol(2_ip)

        else if(words(1)=='PLEPP') then
           if(exists('CONDU')) kfl_plepp_tem = 0
           if(exists('LOWMA')) kfl_plepp_tem = 3
           if(exists('ENTHA')) kfl_plepp_tem = 4

        else if( words(1) == 'NEGAT' ) then
           if(words(2)=='ON   ' ) kfl_negat_tem = 1_ip  ! clip negative overshoots

        else if( words(1) == 'POSIT' ) then
           if(words(2)=='ON   ' ) kfl_posit_tem = 1_ip  ! clip positive overshoots

        end if
     end do
  end if

end subroutine tem_reanut
