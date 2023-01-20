!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_reanut.f90
!> @date    16/11/1966
!> @author  Mariano Vazquez
!> @brief   Read numerical parameters
!> @details Read numerical parameters
!> @}
!------------------------------------------------------------------------
subroutine exm_reanut
  use      def_parame
  use      def_inpout
  use      def_master
  use      def_domain
  use      mod_memchk
  use      def_solver

  use      def_exmedi
  use mod_ecoute, only :  ecoute

  implicit none
  integer(ip) :: istab,iauxi

  if(kfl_paral<=0) then

!
!  Initializations (defaults)
!

     kfl_repro_exm = 0                                ! Res. projection not used 
     kfl_shock_exm = 0                                ! Shock capturing off
     kfl_weigh_exm = 1                                ! dT/dt is in the residual
     kfl_normc_exm = 2                                ! L2 norm for convergence
     kfl_algso_exm = 8                                ! Algebraic solver is GMRES
     kfl_nolim_exm = 0                                ! No non-linear correction method used
     kfl_nolum_exm = 1                                ! Non-linear terms are lumped
     cpu_exmed = 0.0_rp                               ! Total/Solver module CPU

     kfl_timet_exm = 1                                ! Default EXPLICIT time treatment
     kfl_tisch_exm = 1                                ! FFD or BFD
     kfl_tiacc_exm = 1                                ! Default for the explicit case
     kfl_ticel_exm = 1                                ! Forward euler is the default for the cell model 
     kfl_comat_exm = 1                                ! Compute amatr
     
     miinn_exm     = 1                                ! One internal iteration
     msste_exm     = 1                                ! One time substep (i.e., no substepping)
     mnoli_exm     = 1                                ! One non-linear iteration (to compute iap)
     staco_exm(1)  = 4.0_rp                           ! Diffusive term
     staco_exm(2)  = 2.0_rp                           ! Convective term
     staco_exm(3)  = 1.0_rp                           ! Reactive term
     shock_exm     = 0.0_rp                           ! SC parameter
     sstol_exm     = 1.0d-5 !F90                      ! Steady-satate tolerance
     safet_exm     = 1.0_rp                           ! Safety factor
     tnoli_exm     = 0.00001_rp                       ! Tolerance for the non-linear intracellular problem
     cotol_exm     = - 1.0_rp                         ! Subiterations tolerance. Default: NEVER reach convergence
     corat_exm     = - 1.0_rp                         ! Subiterations tolerance ratio. 
     resid_exm     = 0.0_rp
     kfl_adres_exm = 0                                ! Subiterations adaptive ratio flag. DEFAULT: no adaptive. 
     
     dtext_exm     = 0.0_rp                           ! Externally fixed time step
     
     solve_sol => solve(1_ip:)

     !
     ! Reach the section
     !
     rewind(lisda)
     call ecoute('exm_reanut')
     do while(words(1)/='NUMER')
        call ecoute('exm_reanut')
     end do
     !
     ! Begin to read data
     !
     do while(words(1)/='ENDNU')
        call ecoute('exm_reanut')
        if(     words(1)=='STABI') then
           do istab = 1,3
              staco_exm(istab) = param(istab)
           end do
        else if(words(1)=='TIMET') then              ! Explicit (1) or Implicit (2) time treatment
           if(words(2)=='EXPLI') then
              ! Default for the explicit is forward euler
              kfl_timet_exm=1
           else if(words(2)=='IMPLI') then
              ! Default for the implicit is backwards euler
              kfl_timet_exm=2
              if(exists('CRANK')) then
                 kfl_tisch_exm=2
              else if(exists('EULER')) then
                 kfl_tisch_exm=1
              end if
              if(exists('TICEL')) then
                 if(exists('RUNGE')) then
                    kfl_ticel_exm= 2
                    call runend("EXM_REANUT: RUNGE-KUTTA TO BE PROGRAMMED. USE EULER.")
                 else if(exists('EULER')) then
                    kfl_ticel_exm= 1
                 end if

              end if

           end if

        else if(words(1)=='ALGEB') then
           !
           ! Solver
           !
           if (kfl_timet_exm == 1) then              
              kfl_timet_exm = 2   ! this is only meaningful for IMPLICIT, so implicit is forced when ALGEB is read
           end if
           call reasol(1_ip)              

        else if(words(1).eq.'GENER') then
           if(exists('EXPLI')) kfl_genal_exm = 1
           if(exists('DECOU')) kfl_genal_exm = 2
           if(exists('MONOL')) kfl_genal_exm = 3
        else if(words(1).eq.'NONLI') then
           if(exists('LUMPE')) kfl_nolum_exm = 1
        else if(words(1)=='TYPEO') then
           if(words(2)=='RESID') kfl_repro_exm = 1
        else if(words(1)=='SHOCK') then
           if(exists('ISOTR').or.exists('ON   ')) then
              kfl_shock_exm = 1
              !!           shock_exm     = getrea('VALUE',0.0_rp,'#Shock capturing parameter')   ! deprecated
              shock_exm     = getrea('FACTO',0.0_rp,'#Shock capturing parameter')
          else if(exists('ANISO')) then
              kfl_shock_exm = 2
              !!           shock_exm     = getrea('VALUE',0.0_rp,'#Shock capturing parameter')   ! deprecated
              shock_exm     = getrea('FACTO',0.0_rp,'#Shock capturing parameter')
           end if
           if (exists('STATI')) then
              kfl_shock_exm = 11  ! stationary residual, no transient terms
           end if
        else if(words(1)=='TEMPO') then
           if(exists('GALER')) kfl_weigh_exm = 0
!        else if(words(1)=='TIMEA') then
!           kfl_tiacc_exm = int(param(1))
!           kfl_tiacc_exm = 2
        else if(words(1)=='SAFET') then
           safet_exm = param(1)
        else if(words(1)=='SUBST') then  ! substepping
           msste_exm = int(param(1))
        else if(words(1)=='CONST') then  
           kfl_comat_exm = 0                              ! THIS WILL BECOME THE DEFAULT OPTION
           if ( kfl_coupl(ID_SOLIDZ,ID_EXMEDI) >= 1 .or. kfl_coupl(ID_EXMEDI,ID_SOLIDZ) >=1) then              
              if (kfl_gcoup_exm > 0) kfl_comat_exm= 1
           end if
        else if(words(1)=='STEAD') then
           sstol_exm = param(1)
        else if(words(1)=='NORMO') then
           if(exists('L1   ')) then
              kfl_normc_exm = 1
           else if(exists('L-inf')) then
              kfl_normc_exm = 0
           end if
        else if(words(1)=='MAXIM') then
           miinn_exm = int(param(1))
        else if(words(1)=='CONVE') then
           cotol_exm = param(1)

        else if(words(1)=='SUBIT' .or. words(1)=='OUTER') then         !Sub- or outer iterations strategy, NEW WAY!!
           iauxi = 1
           call runend('EXM_REANUT: SUBITERATIONS NOT READY')
           do while (iauxi==1)
              call ecoute('exm_reanut')
              if (words(1)=='ENDSU' .or. words(1)=='ENDOU') iauxi = 0    ! temporary, for back-compatibility
              if (exists('CONVE')) then
                 miinn_exm= getint('ITERA',   1_ip ,'Maximum number of sub-iterations')
                 cotol_exm= getrea('TOLER', -1.0_rp,'Tolerance for the sub-iterations')
                 corat_exm= getrea('RATIO',  0.1_rp,'Tolerance for the sub-iterations')
                 if (exists('ADAPT')) kfl_adres_exm = 1
              end if
           end do

        else if(words(1)=='NLMET') then
           if (exists('ADRIA')) kfl_nolim_exm =  1
        else if(words(1)=='NLITE') then
           mnoli_exm = int(param(1))
        else if(words(1)=='NLTOL') then
           tnoli_exm = param(1)
        else if(words(1)=='ALGEB') then
           if(exists('DIREC')) kfl_algso_exm =  0
           if(exists('CG   ')) kfl_algso_exm =  1
           if(exists('CGNOR')) kfl_algso_exm =  2
           if(exists('BiCG ')) kfl_algso_exm =  3
           if(exists('BiCGW')) kfl_algso_exm =  4
           if(exists('BiCGS')) kfl_algso_exm =  5
           if(exists('TRANS')) kfl_algso_exm =  6
           if(exists('FULLO')) kfl_algso_exm =  7
           if(exists('GMRES')) kfl_algso_exm =  8
           if(exists('FLEXI')) kfl_algso_exm =  9
           if(exists('QUASI')) kfl_algso_exm = 10
        else if(words(1)=='SOLVE') then
           msoit_exm = int(param(1))
        else if(words(1)=='TOLER') then
           solco_exm = param(1)
        else if(words(1)=='KRYLO') then
           nkryd_exm = int(param(1))


        end if
     end do

     if (dtime >= 0._rp ) dtext_exm    = dtime
     if(kfl_genal_exm==1 .or. kfl_genal_exm==2)  then
        nunkn_exm = npoin
        if (kfl_genal_exm==1) kfl_algso_exm =  -2
     else if(kfl_genal_exm==3)  then
        nunkn_exm = ndofn_exm *  npoin
     end if

     if (kfl_timet_exm == 1) then              
        solve_sol(1) % kfl_preco  = SOL_LOCAL_DIAGONAL
     end if


     !if (kfl_cellmod(imate) == 0) kfl_nolum_exm = 0

  end if
  
end subroutine exm_reanut
