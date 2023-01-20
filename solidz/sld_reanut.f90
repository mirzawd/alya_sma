!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_reanut.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   Read numerical treatment
!> @details This routine reads the numerical treatment for solidz module
!> @}
!-----------------------------------------------------------------------

subroutine sld_reanut()
  !
  !.md<module>solidz
  !.md<input>case.sld.dat
  !.md<pos>1
  !.md<sec>
  !
  use def_kintyp,   only : ip, rp
  use def_inpout,   only : words, param, exists, option
  use def_inpout,   only : getint, getrea
  use def_master,   only : solve, solve_sol, INOTSLAVE
  use mod_ecoute,   only : ecoute
  use mod_messages, only : messages_live
  use def_solidz

  implicit none

  external    :: reasol
  external    :: realen
  
  real(rp)    :: vauxi
  integer(ip) :: iauxi

  if( INOTSLAVE ) then
     !
     !  Initializations (defaults)
     !
     kfl_stabi_sld = 0                          ! Stabilization
     kfl_xfeme_sld = 0                          ! XFEM family of enrichment strategies
     kfl_resid_sld = 0                          ! Residual by increment of displacement
     kfl_timet_sld = 1                          ! Time treatment (explicit=1, implicit=2)
     kfl_ninex_sld = 0                          ! Inexact newton (only for implicit cases)
     kfl_tisch_sld = 0                          ! Time integration scheme
     kfl_serei_sld = 0                          ! Selective reduced integration is off
     kfl_limit_sld = 0
     kfl_pseud_sld = 0
     kfl_penal_sld = 0
     kfl_vecto_sld = .false.
     kfl_savdt_sld = 0
     kfl_dttyp_sld = 0
     kfl_celen_sld = 0
     kfl_gdepo     = 1                          ! Compute and store inverse of the deformation gradient
     kfl_ellen_sld = 0                          ! Mininum element lenght by default 
     kfl_safet_table_sld = 0
     minex_sld     = 2                          ! Newton inexact counter
     miinn_sld     = 1                          ! Max # of N-R iterations
     cotol_sld     = 1.0e-5_rp                  ! Convergence tolerance
     safet_sld     = 1.0_rp                     ! Safety factor for time step
     safex_sld     = 1.0_rp                     ! 
     safma_sld     = 1.0e9_rp                   ! Maximum safety factor for time step
     safet_pseud_sld     = 1.0_rp               ! Safety factor for pseudotime step (this value is changed when pseudo is used)
     safet_table_sld = 0.0_rp
     nisaf_sld     = 0
     factor_penal_sld  = 1.0_rp                 ! Penalty factor (this value is changed when penalization is used)
     dafac_sld     = 0.0_rp                     ! Damping factor for automatic stabilization of static problems
     sstol_sld     = 1.0e-5_rp                  ! Steady state tolerance
     masss_sld     = 1.0_rp                     ! Mass scaling factor
     dtmin_sld     = 1.0e-10_rp
     epsex_sld     = 1.0e-8_rp
     !
     ! Reach the section
     !
     call ecoute('sld_reanut')
     do while(words(1)/='NUMER')
        call ecoute('sld_reanut')
     end do
     !
     ! Begin to read data
     !
     !.md# Numerical Treatment Definition
     !.md<code>
     !.md<0>NUMERICAL_TREATMENT
     !
     !.md<com><strong>Numerical treatment</strong>
     !
     do while(words(1)/='ENDNU')
        call ecoute('sld_reanut')

        if(words(1)=='TIMET') then              ! Time treatment
           !
           !.md<1>TIME_TREATMENT:   IMPLICIT | EXPLICIT                      $ Temporal integration scheme
           !.md<field>TIME_TREATMENT
           !.md<com>Set the temporal integration method to IMPLICIT or EXPLICIT methods.
           !.md<com>
           !
           if(words(2)=='EXPLI') then
              kfl_timet_sld=1
           else if(words(2)=='IMPLI') then
              kfl_timet_sld=2
              if(exists('INEXA')) then
                 kfl_ninex_sld= 1
                 epsex_sld= getrea('EPSIL',1.0e-8_rp,'Perturbation to compute the secant stiffness matrix.')
              end if
           end if

        else if(words(1)=='TIMEI') then         ! Time integration scheme
           !
           !.md<1>TIME_INTEGRATION: NEWMARK | TCHAMWA_WIELGOSZ | RUNGE_KUTTA $ Integration scheme method
           !.md<field>TIME_INTEGRATION
           !.md<com>Time integration schemes for the problem.  
           !.md<com>
           !
           if( exists('NEWMA') ) then
              !
              !.md<com>    - NEWMARK: Newmark-Beta Scheme
              !.md<com>The Newmark Beta-method, `NEWMARK, is used by default.`c
              !.md<com>The option `UNDAMPED` uses the undamped trapezoidal rule (`NEWBE=0.25`, `NEWGA=0.5`, `NEWAL=0.0`).
              !.md<com>The option `DAMPED uses the numerically damped integrator with damping proportional to NEWGA-0.5`
              !.md<com>(NEWBE=0.65, NEWGA=0.9, NEWAL=0.0).
              !.md<com>The option `CENTRAL_DIFERENCE uses the explicit central differences method (`NEWBE=0`, `NEWGA=0.5`,`NEWAL=0.0).
              !.md<com>The `CENTRAL_DIFFERENCE` option can only be used for `EXPLICIT` analysis. Moreover, the user can set manually
              !.md<com>`NEWBE=`<tt>real</tt>, `NEWGA=`<tt>real</tt> and `NEWAL=<tt>real</tt>, when one of the previous
              !.md<com>options are missing.
              !.md<com>
              !
              kfl_tisch_sld = 1_ip

              ! Defaults
              if (kfl_timet_sld == SLD_IMPLICIT_SCHEME .and. kfl_timei_sld == SLD_STATIC_PROBLEM ) then      ! Implicit
                 ! This parameters are not used in equilibrium (static) problems
                 tifac_sld(1)  = getrea('NEWBE',0.00_rp,'Beta factor')
                 tifac_sld(2)  = getrea('NEWGA',0.00_rp,'Gamma factor')
                 tifac_sld(3)  = getrea('NEWAL',0.00_rp,'Alpha factor')
              else if (kfl_timet_sld == SLD_IMPLICIT_SCHEME .and. kfl_timei_sld == SLD_DYNAMIC_PROBLEM ) then ! Implicit
                 ! UNDAM
                 tifac_sld(1)  = getrea('NEWBE',0.25_rp,'Beta factor')
                 tifac_sld(2)  = getrea('NEWGA',0.50_rp,'Gamma factor')
                 tifac_sld(3)  = getrea('NEWAL',0.00_rp,'Alpha factor')
                 ! Compatibility check (only when NEWBE is defined by the user)
                 if (tifac_sld(1) == 0.00_rp) call runend("SLD_REANUT: IMPLICIT SCHEMES ARE MEANINGLESS WHEN BETA=0")
              else                                                              ! Explicit
                 ! CDIFF
                 tifac_sld(1)  = getrea('NEWBE',0.00_rp,'Beta factor')
                 tifac_sld(2)  = getrea('NEWGA',0.50_rp,'Gamma factor')
                 tifac_sld(3)  = getrea('NEWAL',0.00_rp,'Alpha factor')
              end if

              if (exists('DAMPE')) then
                 tifac_sld(1)  = 0.65_rp
                 tifac_sld(2)  = 0.90_rp
                 tifac_sld(3)  = 0.00_rp
              else if (exists('UNDAM')) then
                 tifac_sld(1)  = 0.25_rp
                 tifac_sld(2)  = 0.50_rp
                 tifac_sld(3)  = 0.00_rp
              else if (exists('CDIFF')) then
                 tifac_sld(1)  = 0.00_rp
                 tifac_sld(2)  = 0.50_rp
                 tifac_sld(3)  = 0.00_rp
              end if

           else if( exists('TCHAM') ) then
              !
              !.md<com>    - TCHAMWA_WIELGOSZ: Dissipative Tchamwa - Wielgosz explicit scheme
              !.md<com>The `TCHAMWA_WIELGOSZ` is a dissipative explicit scheme. The dissipative parameter 'PHITW' 
              !.md<com>is set to `1.033` by default to quick damping of spurious oscillations associated with instantaneous
              !.md<com>loadings.
              !.md<com>
              !
              kfl_tisch_sld = 2_ip
              tifac_sld(4)  = getrea('PHITW',1.033_rp,'Phi factor')

           else if( exists('RUNGE') ) then
              !
              !.md<com>    - RUNGE: Runge-Kutta Scheme
              !.md<com>The `RUNGE` is a 4th order Runge-Kutta for the resolution of rigid bodies.
              !.md<com>
              !
              kfl_tisch_sld = 3_ip

           else

              call runend('SLD_REANUT: TIME INTEGRATION SCHEME NOT DEFINED')

           end if

        else if( words(1)=='SAFET' ) then         ! Safety factor
           !
           !.md<1>SAFETY_FACTOR= real                                        $ Safety factor
           !.md<field>SAFETY_FACTOR
           !.md<com>Set this parameter equal to the safety factor. This parameter is
           !.md<com>applied to the critical time step. In Explicit analysis is recommended a value smaller
           !.md<com>than one, while in Implicit analysis greater than one. By default, this parameter is set to 1.
           !.md<com>
           !
           safet_sld = param(1)
           !
           safex_sld  = getrea('EXPON',1.0_rp,'#Safety Factor Exponential time function')
           safma_sld  = getrea('MAXIM',1.0e9_rp,'#Maximum safety factor function')
           nisaf_sld  = getint('INITI',0_ip    ,'#Initial step for variable CFL')
           if (exists('TABLE')) then
              iauxi= 0
              call ecoute('sld_reanut')
              do while(words(1)/='ENDTA')
                 iauxi= iauxi+1
                 safet_table_sld(1,iauxi)  = getrea('VALUE',1.0_rp,'#Safety Factor')
                 safet_table_sld(2,iauxi)  = getrea('RESID',0.0_rp,'#Residual')
                 if (iauxi == 10) &
                      call runend('SLD_REANUT: SAFET TABLE MUST BE SMALLER THAN 10.')
                 call ecoute('sld_reanut')
              end do
              safet_sld= safet_table_sld(1,1)  ! starting safety factor value
              kfl_safet_table_sld = iauxi
           end if

        else if( words(1) == 'SAVET' ) then
           !
           !.md<1>SAVE_TIME_STEP:  ON | OFF                                  $ Save time step
           !.md<field>SAVE_TIME_STEP
           !.md<com>Set this parameter to `ON` to save the time step at the begining of the simualtion.
           !.md<com>It avoids the recomputation of the stable time increment for each time step and it may reduce
           !.md<com>the computational cost of the simulation. It can only be used when the 
           !.md<com>strategy for the calculation of the time step uses the underformed shape `REFERENCE`. By default,
           !.md<com>this parameter is deactivated.
           !.md<com>
           !
           if( option('SAVET') ) kfl_savdt_sld = 1
           
        else if( words(1) == 'TIMES' ) then
           !
           !.md<1>TIME_STEP_STRATEGY:  MIN_LENGTH                            $ Time step strategy
           !.md<field>TIME_STEP_STRATEGY:
           !.md<com>Define the time strategy for the problem. By default, the minimum element length method `MIN_LENGHT` is used.
           !.md<com>
           !
           if( words(2) == 'MINLE' ) then
              !
              !.md<com>    - `MIN_LENGHT`: It takes the minimum element length of the model.
              !.md<com>The user can add `REFERENCE` to specify if time step is calculated using the undeformed configuration or 
              !.md<com>`UPDATED` using the deformed configuration.
              !.md<com>
              !.md<com>
              !
              kfl_dttyp_sld = 0
              if(      exists('UPDAT') ) then
                 kfl_celen_sld = 1
              else if( exists('REFER') ) then
                 kfl_celen_sld = 0
              end if
                         
           end if
           
        else if( words(1) == 'DTCRI' ) then
           !
           !.md<1>DTCRI_MIN=  real                                           $ Minimum value for the critical time step  
           !.md<field>DTCRI_MIN:
           !.md<com>Set a mininum critical time step for the simulation to avoid lower values than the computed time step.
           !.md<com>
           !
           dtmin_sld = getrea('DTCRI',1.0e-12_rp,'#Minimum critical time step')

        else if( words(1) == 'ELEME' ) then
           !
           !.md<1>ELEMENT_LENGHT:    MINIMUM | MAXIMUM | VOLUME/AREA | LFACE $ Type of calculation for element characteristic length
           !.md<field>ELEMENT_LENGTH:
           !.md<com>Element characteristic length calculation.
           !.md<com>    - `MINIMUM`: It corresponds to the minimum element edge length.
           !.md<com>    - `MAXIMUM`: It corresponds to the maximum element edge length.
           !.md<com>    - `VOLUME` or `AREA`: It corresponds to the cube or square root for `HEX08` (3-d) and `QUAD4` (2-d) elements
           !.md<com>    respectively.
           !.md<com>    - `LFACE`: It is the volume of the element divided by the largest face area of the element for hexaedrons.
           !.md<com>
           !
           call realen(kfl_ellen_sld)
           
        else if(words(1) == 'STABI') then
           !
           ! ADOC[1]> STABILIZATION: real                         $ Damping factor
           ! ADOC[d]> STABILIZATION: Automatic stabilization with constant damping for static problems.
           !
           kfl_stabi_sld = 1
           dafac_sld = param(1)

        else if(words(1)=='SUBIT') then         !Sub-iterations strategy, NEW WAY!!

           do while (words(1)/='ENDSU')
              call ecoute('sld_reanut')
              if(words(1)  =='PSEUD') then   ! Tau-subiterations (i.e. dual or pseudo time stepping)
                 kfl_pseud_sld = 1
!                 safet_pseud_sld = getrea('SAFET',0.1_rp,'Safety factor for the pseudo-timestep')
                 vauxi = getrea('FACTOR',100.0_rp,'Factor for the pseudo-timestep')
                 safet_pseud_sld = sqrt(1.0_rp / vauxi / vauxi)
              else if(words(1)  =='PENAL') then                        ! Penalization
                 kfl_penal_sld = 1
                 factor_penal_sld = getrea('FACTOR',10.0_rp,'Penalty factor')
              else if(words(1)=='MAXIM') then         ! Linearization iterations
                 miinn_sld  = int(param(1),ip)
              else if(words(1)=='CONVE') then         ! Convergence tolerance
                 cotol_sld = param(1)
              else if(words(1)=='RESID') then         ! Convergence criterion
                 if(exists('DISPL')) then
                    kfl_resid_sld = 0
                 else if(exists('FORCE')) then
                    kfl_resid_sld = 1
                 else if(exists('TOTAL')) then
                    kfl_resid_sld = 2
                 end if
              end if
           end do

        else if(words(1)=='SELEC') then         ! Selective reduced integration enabled
           kfl_serei_sld = 1

        else if( words(1) == 'STEAD' ) then         ! Steady state tolerance
           !
           !.md<1>STEADY_STATE: ON | OFF  TOLERANCE=real                    $ Steady state
           !.md<field>STEADY_STATE
           !.md<com>Steady state in transient problems. The user has to set the tolerance for detection
           !.md<com>of the steady state of the solution. When the time residual calculated between time steps achieves the
           !.md<com>steady state tolerance the analysis stops.
           !.md<com>
           !
           sstol_sld = getrea('STEAD',1.0e-5_rp,'Steady-state tolerance')
           if( exists('OFF  ')) sstol_sld = -1.0_rp
           if( words(2) == 'ON   ' ) then
              if( exists('TOLER')) sstol_sld = getrea('TOLER',1.0e-5_rp,'Steady-state tolerance')
           end if
    
        else if( words(1) == 'NOSTE' ) then         ! Force no steady state check
           !
           !.md<1>NO_STEADY_STATE                                            $ Steady state deactivated
           !.md<field>NO_STEADY_STATE
           !.md<com>Steady state deactivated.
           !.md<com>
           !
           sstol_sld = -1.0_rp

        else if( words(1) == 'MASSS' ) then         ! Mass scaling
           !
           !.md<1>MASS_SCALING= real                                         $ Mass scaling
           !.md<field>MASS_SCALING
           !.md<com>Mass scaling factor for solving quasi-static simulations in Explicit analysis.
           !.md<com>Artificially increaseing the material density, <tt>rho</tt>, by a factor <tt>f^2</tt>
           !.md<com>reduces <tt>n</tt> to <tt>n/f</tt>, just like decreasing <tt>T</tt> to <tt>T/f</tt>.
           !.md<com>This concept reduces the ratio of the event time to time for wave propagation across an
           !.md<com>element while leaving the event time fixed, which allows rate-dependent behavior to be
           !.md<com>included in the analysis. Mass scaling has exactly the same effect on inertia forces as
           !.md<com>speeding up the time of the simulations.
           !.md<com>
           !
           masss_sld = param(1)
           densi_sld(1,1:nmate_sld) = masss_sld*densi_sld(1,1:nmate_sld)

        else if(words(1)=='RESID') then         ! Convergence criterion
           !
           !.md<1>RESIDUAL:         DISPL | FORCE | ENERG | TOTAL            $ Convergence criteria (Implicit analysis)
           !.md<field>RESIDUAL
           !.md<com>This option is used to choose the convergence criteria in the Newton-Raphson iterations.
           !.md<com>It is only applicable to Implicit analysis. By default, the displacement increment error,
           !.md<com>DISPL, is used.
           !.md<com>
           !
           if(exists('DISPL')) then
              !
              !.md<com>    - Set DISPL: Criteria based on the magnitude of the displacement increment.
              !.md<com>
              !
              kfl_resid_sld = 0_ip
           else if(exists('FORCE')) then
              !
              !.md<com>    - Set FORCE: Criteria based on the magnitude of the residual of forces.
              !.md<com>
              !
              kfl_resid_sld = 1_ip
           else if(exists('ENERG')) then
              !
              !.md<com>    - Set ENERGY: Criteria based on an energies.
              !.md<com>
              !
              kfl_resid_sld = 2_ip
           else if(exists('TOTAL') .or. exists('ALL  ')) then
              !
              !.md<com>    - Set ALL: The three criteria have to be fulfilled with the same tolerance.
              !.md<com>
              !
              kfl_resid_sld = 3_ip
           end if

        else if(words(1)=='MAXIM') then         ! Maximum N-R iterations
           !
           !.md<1>MAXIMUM_ITERATION= int                                     $ Max. no. of N-R iterations
           !.md<field>MAXIMUM_ITERAITON
           !.md<com>Maximum number of N-R iterations. This option is only applicable to
           !.md<com>Implicit analysis. By default is set to 1.
           !.md<com>
           !
           miinn_sld = int(param(1),ip)

        else if(words(1)=='CONVE') then         ! Convergence tolerance
           !
           !.md<1>CONVERGENCE_TOL= real                                      $ Convergence tolerance 
           !.md<field>CONVERGENCE_TOL
           !.md<com>Convergence tolerance of the Newton-Raphson iterations. This tolerance is
           !.md<com>the same for all the convergence criteria when TOTAL is used. By default is set to 1e-5.
           !
           cotol_sld = param(1)

        else if(words(1)=='ALGEB') then
           !
           !.md<1>ALGEBRAIC_SOLVER
           !.md<2>...                                                    $ Set solver parameters
           !.md<1>END_ALGEBRAIC_SOLVER
           !.md<field>ALGEBRAIC_SOLVER
           !.md<com>Definition of the algebraic solver.
           !.md<com>
           !
           call reasol(1_ip)

        else if(words(1)=='PRECO') then
           call reasol(2_ip)

        else if(words(1)=='ROTAT')then
           !
           ! ADOC[1]> ROTATION:          On | Off                 $ Rotation to ref. configuration (push forward)
           ! ADOC[d]> ROTATION: Rotation to reference configuration (push forward)
           !
           if(      words(2) == 'ON   ' ) then
              kfl_gdepo = 1_ip
           else if( words(2) == 'OFF  ' ) then
              kfl_gdepo = 0_ip
           end if

        else if(words(1)=='ENRIC') then         ! Enrichement strategy: XFEM/GFEM
           kfl_xfeme_sld= 1
           kfl_xfcra_sld= 1                     ! Local crack definition (as opposed to level set definition)
           if(exists('X-FEM')) kfl_xfeme_sld= 1     ! Default value
           if(exists('LOCAL')) kfl_xfcra_sld= 1     ! Local crack definition
           if(exists('LEVEL')) kfl_xfcra_sld= 2     ! Level set crack definition

        else if(words(1)=='SHELL') then
           do while( words(1) /= 'ENDSH' )
              if(words(1)=='ALGEB') then
                 solve_sol => solve(2:)
                 call reasol(1_ip)
              end if
              call ecoute('SLD_REANUT')
           end do
           solve_sol => solve(1:)

        else if(words(1)=='LIMIT') then
           kfl_limit_sld = 1_ip

        else if(words(1)=='VECTO')then
           if( option('VECTO') ) then
              call messages_live("SOLIDZ: VECTORISED ASSEMBLY")
              kfl_vecto_sld = .true.
           else
              kfl_vecto_sld = .false.
           end if
        end if

     end do
     !
     !.md<0>END_NUMERICAL_TREATMENT
     !.md</code>
     !
     !
     ! Final compatibility corrections
     !
     if (kfl_timet_sld==1 .and. kfl_pseud_sld == 0) then
        miinn_sld     = 1
     end if

  end if

end subroutine sld_reanut
