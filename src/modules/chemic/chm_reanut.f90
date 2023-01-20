!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_reanut
  !-----------------------------------------------------------------------
  !****f* partis/chm_reanut
  ! NAME
  !    chm_reanut
  ! DESCRIPTION
  !    This routine reads the numerical treatment
  ! USES
  !    ecoute
  ! USED BY
  !    chm_turnon
  !***
  !-----------------------------------------------------------------------
  use def_inpout,             only : words, exists, param, getint, getrea
  use def_master,             only : INOTSLAVE, momod, modul, solve
  use def_solver,             only : solve_sol
  use def_chemic,             only : transf_spec_CMC_chm, bemol_chm, chem_int_iZ_CMC_chm, chemical_time_factor, cotol_chm,&
                                     cutof_chm, dampi_chm, dtmax_chm, dtmin_chm, epsht_chm, epstr_chm, kfl_avg_cond_CMC_chm,&
                                     kfl_disable_entpre_chm, kfl_dtcri_chm, kfl_ellen_chm, kfl_entropy_chm, kfl_int_chm,&
                                     kfl_limit_chm, kfl_mesh_interp_CMC_chm, kfl_model_chm, kfl_negat_chm, kfl_norma_chm,&
                                     kfl_normc_chm, kfl_posit_chm, kfl_shock_chm, kfl_solve_cond_CMC_chm, kfl_split_chm,&
                                     kfl_stabi_chm, kfl_taust_chm, kfl_temli_chm, kfl_tiacc_chm, kfl_tibub_chm, kfl_timei_chm,&
                                     kfl_tisch_chm, kfl_transfer_condField_CMC_chm, kfl_wallc_chm, kfl_warni_chm, nclas_chm,&
                                     nspec_transf_CMC_chm, nZ_chm_int_CMC_chm, nZ_CMC_chm, order_RK_CMC_chm, relax_chm, safet_chm,&
                                     shock_chm, sstol_chm, strec_chm, temli_chm, transf_entha_CMC_chm, Zs_CMC_chm, staco_chm,&
                                     extr_Z_chem_integr_CMC_chm, posZ_chem_integr_CMC_chm
  use def_kintyp,             only : ip, rp
  use mod_ecoute,             only : ecoute
  use mod_chm_operations_CMC, only : chm_initial_actions_reanut_CMC, &
                                     chm_chemistry_limits_CMC, &
                                     chm_construct_vector_chem_integ_CMC
  implicit none
  integer(ip) :: ivari, ii, jj, aux_nZ, kfl_chm_int = 0_ip
  integer(ip), pointer ::  aux_Z(:)

  external             :: reatau
  external             :: realen
  external             :: reasol
  external             :: runend
  external             :: chm_memnut

  if( INOTSLAVE ) then
     !
     !  Initializations (defaults)
     !
     kfl_ellen_chm = 0                                ! Minimum element length
     kfl_taust_chm = 1                                ! Tau strategy
     kfl_shock_chm = 0                                ! Shock capturing off
     kfl_stabi_chm = 0                                ! Galerkin
     kfl_limit_chm = 0                                ! No limiter
     kfl_wallc_chm = 0                                ! Flamelet combustion model correction source term at walls (=0 OFF, =1 ON)
     kfl_tibub_chm = 0                                ! Time integration of Bubble function
     kfl_tiacc_chm = 1                                ! First order time integ.
     kfl_tisch_chm = 1                                ! Trapezoidal rule
     kfl_normc_chm = 2                                ! L2 norm for convergence
     kfl_dtcri_chm = 1                                ! Time step criteria
     kfl_negat_chm = 0                                ! Strategy for negative concentrations
     kfl_posit_chm = 0                                ! Strategy for too positive concentrations
     kfl_warni_chm = 1                                ! DEfault warn about points where mass sums to zero
     kfl_split_chm = -1                               ! Splitting algrithm, ORDER = 1,2
     kfl_int_chm = -1                                 ! Chemical Source integration (Default CVODE)
     kfl_transfer_condField_CMC_chm = 0_ip            ! 0: do not transfer conditional fields, 1: transfer conditional fields
     staco_chm(1)  = 1.0_rp                           ! Diffusive term
     staco_chm(2)  = 1.0_rp                           ! Convective term
     staco_chm(3)  = 1.0_rp                           ! Reactive term
     shock_chm     = 0.0_rp                           ! SC parameter
     bemol_chm     = 0.0_rp                           ! Bemol (convectice term)
     temli_chm     = 0.0_rp                           ! T limiter to compute reaction rates
     cotol_chm     = 1.0e-3_rp                        ! Convergence tolerance
     safet_chm     = 1.0_rp                           ! Safety factor
     chemical_time_factor = 1.0_rp                    ! Safety factor exclusively for the source term
     cutof_chm     = 1.0e-8_rp                        ! Concentration cutoff for critical time computation
     sstol_chm     = -1.0e-5_rp                       ! Steady-state tolerance
     strec_chm     = 2.0_rp                           ! Adaptive dt: Streatching factor
     dampi_chm     = 2.0_rp                           ! Adaptive dt: damping
     epsht_chm     = 0.025_rp                         ! Adaptive dt: eps_R
     epstr_chm     = 0.025_rp                         ! Adaptive dt: eps_A
     dtmin_chm     = 1.0e-12_rp                       ! Minimum time step
     dtmax_chm     = 1.0e12_rp                        ! Maximum time step
     kfl_temli_chm = 0                                ! Flag to activate a T limiter to compute reaction rates
     relax_chm = 1.0_rp                               ! Relaxation

     kfl_entropy_chm = 0                              ! entropy viscosity stablization method
     kfl_disable_entpre_chm = 0                       ! Option to disable prediction step in entropy viscosity model

     kfl_mesh_interp_CMC_chm = 0_ip                   ! 0: same meshes for CMC and CFD, 1: different meshes
     kfl_avg_cond_CMC_chm = 1_ip                      ! 0: averages with unconditional fields, 1: averages with conditional fields
     extr_Z_chem_integr_CMC_chm(1) = 0.0_rp
     extr_Z_chem_integr_CMC_chm(2) = 1.0_rp
     nspec_transf_CMC_chm = 0_ip
     transf_entha_CMC_chm = 0_ip

     solve_sol     => solve                           ! Solver type
     ivari         =  1
     !
     ! Reach the section
     !
     call ecoute('chm_reanut')
     do while(words(1)/='NUMER')
        call ecoute('chm_reanut')
     end do
     !
     ! Begin to read data
     !
     do while(words(1)/='ENDNU')
        call ecoute('chm_reanut')

        !----------------------------------------------------------------
        !
        ! Stabilization strategy
        !
        !----------------------------------------------------------------

        if( words(1) == 'TAUST' ) then
           call reatau(kfl_taust_chm)

        else if( words(1) == 'STRAT' .or. words(1) == 'STABI' ) then
           if( words(2) == 'GALER' .or. words(2) == 'OFF  ' ) then
              kfl_stabi_chm = -2
           else if ( words(2) == 'ENTRO') then
              if (words(3) == 'ENTPR') then
                 if (words(4) == 'OFF') then
                    kfl_disable_entpre_chm = 1
                 else if (words(4) == 'ON') then
                    kfl_disable_entpre_chm = 0
                 end if
              end if
              kfl_stabi_chm = -2
              kfl_entropy_chm = 1
           else if( words(2) == 'SU   ' .or. words(2) == 'FIRST' ) then
              kfl_stabi_chm = -1
           else if( words(2) == 'ASGS ' ) then
              kfl_stabi_chm =  0
           else if( words(2) == 'FULLO' ) then
              kfl_stabi_chm =  1
           else if( words(2) == 'OSS  ' .or. words(2) == 'AOSS ') then
              kfl_stabi_chm =  2
              if( exists('NOLIM') ) kfl_limit_chm =  0  ! No limiter
              if( exists('SOTO ') ) kfl_limit_chm =  1  ! Soto limiter
              if( exists('DIFFU') ) kfl_limit_chm =  2  ! Very diffusive limiter
              if( exists('FIRST') ) kfl_limit_chm = -1  ! First order
           else if( words(2) == 'AROSS' ) then
              kfl_stabi_chm =  3
              if( exists('NOLIM') ) kfl_limit_chm =  0  ! No limiter
              if( exists('SOTO ') ) kfl_limit_chm =  1  ! Soto limiter
              if( exists('DIFFU') ) kfl_limit_chm =  2  ! Very diffusive limiter
              if( exists('FIRST') ) kfl_limit_chm = -1  ! First order
           else if( words(2) == 'MARGA' ) then
              kfl_stabi_chm =  4
           else if( words(2) == 'SUPG ' ) then
              kfl_stabi_chm =  -3
           else if( words(2) == 'BUBBL' ) then
              kfl_stabi_chm =  -4
              if( words(3) == 'TIMET' ) kfl_tibub_chm = 1
           else if(words(2)=='CONST' ) then ! reads stab constants
              staco_chm(1) = param(2)
              staco_chm(2) = param(3)
              staco_chm(3) = param(4)
           end if

        else if( words(1) == 'ELEME' ) then
           call realen(kfl_ellen_chm)

        else if( words(1) == 'SPLIT' ) then
              kfl_split_chm = getint('ORDER',kfl_split_chm,'#Splitting algorithm order')
              write(momod(modul) % lun_outpu,*)''
              write(momod(modul) % lun_outpu,*)'SPLITTING ORDER = ',kfl_split_chm
              write(momod(modul) % lun_outpu,*)''

        else if( words(1) == 'CMCMO' ) then
           !
           ! Data for CMC model
           !
           call ecoute('chm_reanut')

           CMC_model: do while(words(1)/='ENDCM')

              if( words(1) == 'MFSPL' ) then  ! RK order for mixture fraction splitting
                 order_RK_CMC_chm = getint('MFSPL',3_ip,'#Order for the Runge-Kutta applied to mixture fraction diffusion')
                 write(momod(modul) % lun_outpu,*)''
                 write(momod(modul) % lun_outpu,*)'SPLITTING ORDER FOR MIXTURE FRACTION = ', order_RK_CMC_chm
                 write(momod(modul) % lun_outpu,*)''

                 if (order_RK_CMC_chm <= 0_ip .or. order_RK_CMC_chm >= 5_ip)  call runend('CHEMIC REANUT: order not available for&
                    & RK in mixt. frac. diffusion')

              else if( words(1) == 'MESHI' ) then  ! Mesh interpolation in CMC
                 if (words(2) == 'YES  ')   kfl_mesh_interp_CMC_chm = 1_ip
                 if (kfl_mesh_interp_CMC_chm == 1_ip) then
                    if (words(3) == 'UNCON')  kfl_avg_cond_CMC_chm = 0_ip
                 end if

              else if( words(1) == 'CHEMI' ) then
                 call ecoute('chm_reanut')
                 do while (words(1) /= 'ENDCH' )
                    if (words(1) == 'MINIM') then
                       extr_Z_chem_integr_CMC_chm(1) = param(1)
                       extr_Z_chem_integr_CMC_chm(1) = max(0.0_rp, extr_Z_chem_integr_CMC_chm(1))

                    else if (words(1) == 'MAXIM') then
                       extr_Z_chem_integr_CMC_chm(2) = param(1)
                       extr_Z_chem_integr_CMC_chm(2) = min(Zs_CMC_chm, extr_Z_chem_integr_CMC_chm(2))

                    else if (words(1) == 'MIXTU') then
                       if ( words(2) == 'TOTAL' ) then
                           aux_nZ = getint('TOTAL',0_ip,'#Number of mixture fractions to be chemically integrated')
                           if (aux_nZ >= nZ_CMC_chm)  call runend('CHEMIC REANUT: number of mixture fractions to be chemically&
                               & integrated is equal or exceeds the number of total mixture fractions')
                       else
                           call runend('CHEMIC REAOUS: number of mixture fractions to be chemically integrated not given')
                       end if
                       if ( words(3) == 'POSIT' .or. words(3) == 'NOPOS' ) then
                          if ( words(3) == 'POSIT')   kfl_chm_int = 1_ip
                          allocate(aux_Z(aux_nZ))
                          do ii = 1, aux_nZ
                             if (param(ii+2)<1 .or. param(ii+2)>real(nZ_CMC_chm,rp)) then
                                call runend('CHEMIC REAOUS: position of mixture fractions to be post-processed not valid')
                             else
                                aux_Z(ii) = int(param(ii+2_ip), kind=ip)
                             end if
                          end do
                       else
                          call runend('CHEMIC REAOUS: position of mixture fractions to be chemically integrated not given')
                       end if

                    end if

                    call ecoute('chm_reanut')
                 end do

              else if ( words(1) == 'TRANS' ) then
                 kfl_transfer_condField_CMC_chm = 1_ip
                 call ecoute('chm_reanut')
                 do while (words(1) /= 'ENDTR')
                    if ( words(1) == 'SPECI' ) then
                       if ( words(2) == 'TOTAL' ) then
                          nspec_transf_CMC_chm = getint('TOTAL',0_ip,'#Number of species to be transferred from CFD to CMC')
                          if (nspec_transf_CMC_chm < 0_ip .or. nspec_transf_CMC_chm >= nclas_chm) &
                             call runend('CHEMIC REANUT: number of species to be transferred not valid')
                          call chm_memnut(2_ip)
                          if ( words(3) == 'POSIT' ) then
                             do ii = 1, nspec_transf_CMC_chm
                                if (param(ii+2)<real(1,rp) .or. param(ii+2)>real(nclas_chm,rp)) then
                                   call runend('CHEMIC REAOUS: number of species to be post-processed not valid')
                                else
                                   transf_spec_CMC_chm(ii) = int(param(ii+2_ip), kind=ip)
                                end if
                             end do
                          else
                             call runend('CHEMIC REANUT: position of species to be transferred not given')
                          end if
                       else
                          nspec_transf_CMC_chm = nclas_chm
                          call chm_memnut(2_ip)
                          transf_spec_CMC_chm = (/( ii, ii=1_ip, nclas_chm )/)
                       end if
                    else if ( words(1) == 'ENTHA' ) then
                       transf_entha_CMC_chm = 1_ip
                    end if
                    call ecoute('chm_reanut')
                 end do
              end if

              call ecoute('chm_reanut')

           end do CMC_model


        else if( words(1) == 'NEGAT' ) then
           if(words(2)=='ON   ' ) then
              if (exists('LAGRA')) then
                 kfl_negat_chm = 0_ip  ! no negative strategy
                 kfl_norma_chm = -2_ip  ! Lagrangian multiplier normalization
              else
                 kfl_negat_chm = 1_ip  ! take previous steps
              endif
           else
              kfl_negat_chm = 0_ip     ! Leave as is, no strategy
           endif

        else if( words(1) == 'POSIT' ) then
           if(words(2)=='ON   ' ) then
              kfl_posit_chm = 1
           end if

        else if( words(1) == 'WARNI' ) then
           if(words(2)=='ON   ' ) then
              kfl_warni_chm = 1
           else if(words(2)=='OFF  ' ) then
              kfl_warni_chm = 0
           end if

        else if( words(1) == 'RELAX' ) then
           if(words(2)=='ON   ' ) then
              relax_chm = getrea('PARAM',0.5_rp,     '#Relaxation factor for update')
           else
              relax_chm = 1.0_rp
           endif

        else if( words(1) == 'ELEME' ) then
           call realen(kfl_ellen_chm)
           if(words(2)=='NEW  ') kfl_ellen_chm=-1

        else if( words(1) == 'TIMEI' ) then
           if(exists('TRAPE').or. exists('BDF  ') ) call runend('In chm_reanut: implicit schemes not available anymore')

           if(exists('ADAMS')) then
              kfl_tisch_chm    = 3
           end if

           if(exists('RUNGE')) then
              kfl_tisch_chm    = 4
           end if

           kfl_tiacc_chm = getint('ORDER',1_ip,'#Time integration order')

        else if( words(1) == 'CHEMI' ) then
           if(exists('CVODE')) then
              kfl_int_chm = -1_ip
              write(momod(modul) % lun_outpu,*)''
              write(momod(modul) % lun_outpu,*)'CVODE Default integration'
              write(momod(modul) % lun_outpu,*)''
           else if(exists('ODEPI')) then
              kfl_int_chm = 1_ip
              write(momod(modul) % lun_outpu,*)''
              write(momod(modul) % lun_outpu,*)'ODEPIM Integration'
              write(momod(modul) % lun_outpu,*)''
           else if(exists('PYJAC')) then
              kfl_int_chm = 2_ip
              write(momod(modul) % lun_outpu,*)''
              write(momod(modul) % lun_outpu,*)'CVODE + PYJAC Integration'
              write(momod(modul) % lun_outpu,*)''
           else if(exists('ROSEN')) then
              kfl_int_chm = 3_ip
              write(momod(modul) % lun_outpu,*)''
              write(momod(modul) % lun_outpu,*)'Rosenbrock4 + PYJAC Integration'
              write(momod(modul) % lun_outpu,*)''
           end if

        else if( words(1) == 'TIMEA' ) then
           kfl_tiacc_chm = int(param(1),ip)
           if(kfl_timei_chm==0) kfl_tiacc_chm = 1

        else if( words(1) == 'CONVE' ) then
           cotol_chm = getrea('CONVE',1.0e-3_rp,'#CONVERGENCE TOLERANCE')

        else if( words(1) == 'SAFET' ) then
           safet_chm = param(1)
           if (exists('SOURC')) chemical_time_factor = getrea('SOURC',1.0_rp,'#Source term safety factor')

        else if( words(1) == 'CONCE' ) then
           cutof_chm = param(1)

        else if( words(1) == 'STEAD' ) then
           sstol_chm = param(1)

        else if (words(1) == 'WALLC' ) then
            if(words(2)=='ON   ' )  kfl_wallc_chm = 1_ip

        else if( words(1) == 'NORMO' ) then
           if(exists('L1   ')) then
              kfl_normc_chm = 1
           else if(exists('L-INF')) then
              kfl_normc_chm = 0
           else if(exists('LINF ')) then
              kfl_normc_chm = 0
           else if(exists('L2   ')) then
              kfl_normc_chm = 2
           else if(exists('ALGEB')) then
              kfl_normc_chm = 3
           end if

        else if( words(1) == 'BEMOL' ) then
           bemol_chm = getrea('BEMOL',0.0_rp,'#Bemol of convective term')

        else if( words(1) == 'CONSI' ) then
           !
           ! Consistent matrix
           !
           ivari = 2

        else if( words(1) == 'ALGEB' ) then
           solve_sol => solve(ivari:)
           call reasol(1_ip)

        else if( words(1) == 'LIMIT' ) then
           kfl_temli_chm = 1
           temli_chm = getrea('TEMPE',-1000.0_rp,'#Temperature limiter to compute reaction rates')

        else if( words(1) == 'PRECO' ) then
           call reasol(2_ip)

        end if
     end do


  if (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1) then
     call chm_chemistry_limits_CMC

     if (.not. associated(aux_Z)) then
        nZ_chm_int_CMC_chm = posZ_chem_integr_CMC_chm(2) - posZ_chem_integr_CMC_chm(1) + 1_ip
        call chm_memnut(1_ip)
        ! Assign
        chem_int_iZ_CMC_chm = (/( jj, jj=posZ_chem_integr_CMC_chm(1), posZ_chem_integr_CMC_chm(2) )/)
     else
        call chm_construct_vector_chem_integ_CMC(aux_nZ, aux_Z, kfl_chm_int)
     end if

     call chm_initial_actions_reanut_CMC
  end if

  end if

end subroutine chm_reanut
